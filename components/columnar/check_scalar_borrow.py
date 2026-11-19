#!/usr/bin/env python3
"""Compare actual borrowed linked columns with the original Awkward projection."""
from collections import Counter
import json
from pathlib import Path
import sys
from unittest.mock import patch

import awkward as ak
import numpy as np
import typer
import uproot

from atlas_met_columnar import ASSOC, MemberwiseJets, Native, pack_raw, read_associations
from atlas_met_prepare import Data, prepare_chunk
from met_demo import UUID, projected_fields, read_branch_names, read_owner_keys


def main(input_file: Path, runtime: Path, output: Path, entries: int = typer.Option(2000,min=1)):
    native, payloads = runtime/"native",runtime/"payloads"
    keys = read_owner_keys(runtime)
    bridge = Native(native/"libatlas_met_columnar.so")
    context = bridge.open(native/"libatlas_met_maker_provider.so",payloads/"maker.json",payloads/"nnjvt.json")
    try:
        bridge.configure(context,payloads/"soft.json",payloads/"significance.json",
                         native/"libmet_preparation_provider.so",payloads/"preparation.json")
        with uproot.open(input_file,handler=uproot.source.file.MemmapSource,
                         num_workers=1,num_fallback_workers=1) as original:
            if str(original.file.uuid) != UUID:
                raise ValueError("Unexpected original file UUID")
            tree = original["CollectionTree"]
            if entries > tree.num_entries:
                raise ValueError("Requested entries exceed the original tree")
            fields = projected_fields(read_branch_names(tree))
            raw = tree.arrays(fields,entry_start=0,entry_stop=entries,library="ak",how=dict)
            associations = read_associations(tree,0,entries)
            jet_links = MemberwiseJets(bridge,tree[ASSOC][ASSOC+"jetLink"]).read(0,entries)
            candidate = pack_raw(raw,associations,jet_links)
            generic = dict(candidate)

            borrowed = {}
            hits = Counter()
            misses = Counter()
            cache_size = 0
            access = Data.scalar_view

            def observe(data,name,target,rows):
                nonlocal cache_size
                view = access(data,name,target,rows)
                cache_size = max(cache_size,sys.getsizeof(data.scalar_views))
                if view is None:
                    misses[name] += 1
                else:
                    if view is not ak.to_layout(data.raw[name]).content.data:
                        raise AssertionError(f"Cached array was copied: {name}")
                    borrowed[name] = {"nbytes":view.nbytes,"dtype":view.dtype.str,"elements":view.size}
                    hits[name] += 1
                    hits[name+":rows"] += rows.size
                return view

            evaluate = lambda name, values: bridge.prepare(context,name,values)
            with patch.object(Data,"scalar_view",observe):
                prepare_chunk(candidate,raw,keys,evaluate,capture_preparation=True)
            with patch.object(Data,"scalar_view",return_value=None):
                prepare_chunk(generic,raw,keys,evaluate,capture_preparation=True)

        if candidate.keys() != generic.keys():
            raise AssertionError("Projection branch inventory changed")
        for name,values in candidate.items():
            stock = generic[name]
            if values.shape != stock.shape or values.dtype != stock.dtype:
                raise AssertionError(f"Projection shape or dtype changed: {name}")
            if values.dtype.hasobject:
                np.testing.assert_array_equal(values,stock)
            elif values.tobytes() != stock.tobytes():
                raise AssertionError(f"Projection bits changed: {name}")
        selected = {role:int(np.count_nonzero(candidate[role+".selected"]))
                    for role in ("electron","photon","muon","tau")}
        report = {
            "entries":entries,"input_uuid":UUID,"exact_prepared_transport_columns":len(candidate),
            "selected_by_species":selected,
            "fast_projection_calls":sum(n for name,n in hits.items() if not name.endswith(":rows")),
            "fast_gathered_rows":sum(n for name,n in hits.items() if name.endswith(":rows")),
            "fast_calls_and_rows_by_leaf":dict(hits),"attempted_scalar_fallbacks_by_leaf":dict(misses),
            "borrowed_leaves":borrowed,"borrowed_payload_bytes":sum(item["nbytes"] for item in borrowed.values()),
            "additional_cached_payload_bytes":0,"shallow_cached_dictionary_bytes":cache_size,
        }
        with output.open("x") as handle:
            json.dump(report,handle,indent=2)
            handle.write("\n")
        if not hits or not all(selected.values()):
            raise ValueError("Raw witness lacks verified borrows or a selected species")
    finally:
        bridge.lib.met_columnar_close(context)


if __name__ == "__main__":
    typer.run(main)
