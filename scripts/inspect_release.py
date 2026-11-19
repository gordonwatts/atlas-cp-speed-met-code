"""Observe the original xAOD containers and default release MET systematics."""
from __future__ import annotations

import json
import os
from pathlib import Path
import sys

import typer

app = typer.Typer(add_completion=False)


def inspect_link(link) -> dict:
    valid = bool(link.isValid())
    return {
        "persistent_key": int(link.persKey()),
        "persistent_index": int(link.persIndex()),
        "key": int(link.key()),
        "index": int(link.index()),
        "valid": valid,
        "container": str(link.dataID()) if valid else None,
    }


@app.command()
def inspect(
    input_file: Path,
    output_file: Path,
    events: int = typer.Option(16, min=1, max=20000),
) -> None:
    import ROOT

    ROOT.xAOD.Init().ignore()
    source = ROOT.TFile.Open(str(input_file), "READ")
    if not source or source.IsZombie():
        raise RuntimeError(f"Cannot open {input_file}")
    original = source.Get("CollectionTree")
    identity = {
        "path": str(input_file),
        "resolved_path": str(input_file.resolve()),
        "uuid": str(source.GetUUID().AsString()),
        "entries": int(original.GetEntries()),
        "bytes": int(source.GetSize()),
        "filesystem_bytes": input_file.stat().st_size,
    }
    tool = ROOT.met.METSystematicsTool("IntakeMETSystematicsTool")
    if not tool.initialize().isSuccess():
        raise RuntimeError("Default release METSystematicsTool did not initialize")
    report = {
        "input": identity,
        "release": os.environ["AtlasVersion"],
        "python": sys.version,
        "root": str(ROOT.gROOT.GetVersion()),
        "root_implicit_mt": bool(ROOT.IsImplicitMTEnabled()),
        "cpu_affinity": sorted(os.sched_getaffinity(0)),
        "process_threads": len(list(Path("/proc/self/task").iterdir())),
        "recommended": sorted(str(state.name()) for state in tool.recommendedSystematics()),
        "affecting": sorted(str(state.name()) for state in tool.affectingSystematics()),
        "events": [],
    }
    transient = ROOT.xAOD.MakeTransientTree(source, "CollectionTree", ROOT.xAOD.TEvent.kClassAccess)
    for entry in range(events):
        if transient.GetEntry(entry) < 0:
            raise RuntimeError(f"Failed to read entry {entry}")
        row = {
            "entry": entry,
            "event_number": int(transient.EventInfo.eventNumber()),
            "jet_count": len(transient.AnalysisJets),
            "associations": [
                {
                    "is_misc": bool(assoc.isMisc()),
                    "jet": inspect_link(assoc.jetLink()),
                    "objects": [inspect_link(link) for link in assoc.objectLinks()],
                }
                for assoc in transient.METAssoc_AnalysisMET
            ],
            "core_terms": [str(term.name()) for term in transient.MET_Core_AnalysisMET],
            "truth_terms": [str(term.name()) for term in transient.MET_Truth],
        }
        report["events"].append(row)
    with output_file.open("x") as output:
        json.dump(report, output, indent=2, allow_nan=False)
        output.write("\n")
    ROOT.xAOD.ClearTransientTrees()
    source.Close()


if __name__ == "__main__":
    app()
