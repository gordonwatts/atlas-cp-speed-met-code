#!/usr/bin/env python3
"""Create the pinned CorrectionLib event operation from the lowered release cut map."""
import json
from pathlib import Path

import typer


def main(lowered: Path, output: Path) -> None:
    data = json.loads(lowered.read_text(encoding="utf-8"))
    if data["recipe"] != "AnalysisBase-25.2.117-PFlow-Tight-PVSoftTrk" or data["working_point"] != "FixedEffPt":
        raise ValueError("Wrong recipe or working point")
    pt = data["pt_edges_mev"]
    eta = data["eta_edges"]
    packed = " ".join(map(str, [len(pt), len(eta), *pt, *eta, *data["cuts"]]))
    result = {
        "schema_version": 3,
        "corrections": [{
            "name": "build_met_terms", "version": 1,
            "inputs": [{"name": "invocation", "type": "int"}],
            "output": {"name": "bound_result_written", "type": "real"},
            "data": {"nodetype": "extension", "type": "atlas.met.build_met_terms", "version": 1,
                     "config": {"recipe": data["recipe"], "source_sha256": data["source_sha256"],
                                "cuts": packed, "min_pt_mev": 20000.0, "max_pt_mev": 60000.0,
                                "min_abs_eta": -1.0, "max_abs_eta": 2.5}},
        }],
    }
    with output.open("x", encoding="utf-8") as handle:
        json.dump(result, handle, allow_nan=False)
        handle.write("\n")


if __name__ == "__main__":
    typer.run(main)
