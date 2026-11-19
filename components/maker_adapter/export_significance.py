#!/usr/bin/env python3
"""Bind the four peer-exported tables to the MET significance event node."""
import json
from pathlib import Path

import typer


def main(exported_primitives: Path, output: Path) -> None:
    data = json.loads(exported_primitives.read_text(encoding="utf-8"))
    tables = {}
    for kind in ("jet", "egamma", "muon", "tau"):
        op = next(item for item in data["corrections"] if item["name"] == f"{kind}_resolution")
        table = Path(op["data"]["config"]["table"])
        if not table.is_absolute():
            table = exported_primitives.resolve().parent / table
        tables[kind] = str(table.resolve(strict=True))
    payload = {
        "schema_version": 3,
        "corrections": [{
            "name": "met_significance_event", "version": 1,
            "inputs": [{"name": "invocation", "type": "int"}],
            "output": {"name": "bound_result_written", "type": "real"},
            "data": {"nodetype": "extension", "type": "atlas.met.significance_event",
                     "version": 1, "config": tables},
        }],
    }
    with output.open("x", encoding="utf-8") as handle:
        json.dump(payload, handle, allow_nan=False)
        handle.write("\n")


if __name__ == "__main__":
    typer.run(main)
