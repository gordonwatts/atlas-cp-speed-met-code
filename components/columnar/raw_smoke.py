"""Real PHYSLITE original-owner/association/JVT transport smoke."""
from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np
import typer
import uproot

from atlas_met_columnar import ASSOC, RAW_BRANCHES, MemberwiseJets, Native, pack_raw, read_associations

UUID = "79ee6652-311f-11f1-aa37-0445010abeef"


def main(input_file: Path, output: Path, native: Path, provider: Path, maker: Path, nnjvt: Path,
         entries: int = typer.Option(16, min=1, max=20000),
         chunk_size: int = typer.Option(16, min=1)) -> None:
    if sys.version_info[:2] != (3, 14):
        raise ValueError("This standalone diagnostic requires CPython 3.14")
    bridge = Native(native)
    context = bridge.open(provider, maker, nnjvt)
    terms, flags, scores, decisions, keys, indices = [], [], [], [], [], []
    try:
        with uproot.open(input_file, handler=uproot.source.file.MemmapSource,
                         num_workers=1, num_fallback_workers=1) as root:
            if str(root.file.uuid) != UUID:
                raise ValueError("Original association owner keys are pinned to the campaign input")
            tree = root["CollectionTree"]
            if entries > tree.num_entries:
                raise ValueError("Input has fewer entries than requested")
            link_reader = MemberwiseJets(bridge, tree[ASSOC][ASSOC+"jetLink"])
            for start in range(0, entries, chunk_size):
                stop = min(start+chunk_size, entries)
                direct = tree.arrays(RAW_BRANCHES, entry_start=start, entry_stop=stop, library="ak", how=dict)
                association = read_associations(tree, start, stop)
                memberwise = link_reader.read(start, stop)
                buffers = pack_raw(direct, association, memberwise)
                values = bridge.raw_maker(context, buffers, stop-start)
                for column, value in zip((terms, flags, scores, decisions), values, strict=True):
                    column.append(value)
                keys.append(memberwise[1])
                indices.append(memberwise[2])
    finally:
        bridge.lib.met_columnar_close(context)
    with output.open("xb") as handle:
        np.savez(handle, term_vectors=np.concatenate(terms), flag_masks=np.concatenate(flags),
                 nnjvt=np.concatenate(scores), jet_jvt=np.concatenate(decisions),
                 association_jet_key=np.concatenate(keys), association_jet_index=np.concatenate(indices))
    print(json.dumps({"stage": "original-association-nnjvt-empty-hard-views", "entries": entries,
                      "uuid": UUID, "output": str(output)}))


if __name__ == "__main__":
    typer.run(main)
