# Copyright (C) 2002-2026 CERN for the benefit of the ATLAS collaboration
# SPDX-License-Identifier: Apache-2.0
"""Export the resolved AnalysisBase 25.2.117 MET NNJVT cut data."""

import hashlib
import json
import math
from pathlib import Path
import re
import struct
from typing import Annotated

import typer


def float32(value: int | float) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"expected a numeric JSON value, got {value!r}")
    rounded = struct.unpack("<f", struct.pack("<f", float(value)))[0]
    if not math.isfinite(rounded):
        raise ValueError(f"nonfinite float32 {value}")
    return rounded


def payload(source: bytes, source_name: str) -> dict:
    original = json.loads(source)
    pt = [float32(float32(x) * 1000.0) for x in original["ptbin_edges"]]
    eta = [float32(x) for x in original["etabin_edges"]]
    for name, edges in (("pt", pt), ("eta", eta)):
        if len(edges) < 2 or any(a >= b for a, b in zip(edges, edges[1:])):
            raise ValueError(f"invalid {name} edges after source float32 conversion")
    num_pt = len(pt) - 1
    num_eta = len(eta) - 1
    cuts = [0.0] * (num_pt * num_eta)
    for key, value in sorted(original["cuts"].items()):
        match = re.fullmatch(r"\((\d+),\s*(\d+)\)", key, re.ASCII)
        if match is None:
            raise ValueError(f"invalid cut bin {key!r}")
        i, j = (int(item) for item in match.groups())
        if i >= num_pt or j >= num_eta:
            raise ValueError(f"out-of-range cut bin {key!r}")
        cuts[i * num_eta + j] = float32(value)
    return {
        "recipe": "AnalysisBase-25.2.117-PFlow-Tight-PVSoftTrk",
        "working_point": "FixedEffPt",
        "source_file": source_name,
        "source_sha256": hashlib.sha256(source).hexdigest(),
        "pt_edges_mev": pt,
        "eta_edges": eta,
        "cuts": cuts,
    }


def main(
    source: Annotated[Path, typer.Argument(exists=True, dir_okay=False, readable=True)],
    output: Annotated[Path, typer.Argument(dir_okay=False)],
) -> None:
    """Export the resolved release JSON to the common provider's typed table."""
    raw = source.read_bytes()
    table = payload(raw, str(source.resolve()))
    with output.open("x", encoding="utf-8") as stream:
        json.dump(table, stream, indent=2, allow_nan=False)
        stream.write("\n")


if __name__ == "__main__":
    typer.run(main)
