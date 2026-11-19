#!/usr/bin/env python3
"""Check to which library each live execute symbol was bound and which ran."""
from collections import Counter
import json
from pathlib import Path
import re

import typer


def main(run: Path, log: Path, expected_calls: int | None = None) -> None:
    binding = json.loads((run / "binding-map.json").read_text(encoding="utf-8"))
    stage = binding["stage"]
    if stage == "maker_only":
        symbols, operations = {"maker"}, {"maker", "stock_soft_diagnostic"}
    elif stage == "full":
        symbols, operations = {"maker", "final", "significance"}, {"maker", "soft", "final", "significance"}
    else:
        raise ValueError(f"unexpected framework stage {stage}")
    if set(binding["symbols"]) != symbols:
        raise ValueError(f"wrong live execute symbols for {stage}")
    if any(row["owner"] != binding["bridge"] for row in binding["symbols"].values()):
        raise ValueError("incorrect live execute owner")
    if binding["worker_threads"] != 1 or len(binding["cpu_affinity"]) != 1:
        raise ValueError("the audited process does not have a single thread and CPU")
    counts = Counter()
    for operation, value in re.findall(r"^MET_B_DISPATCH (\S+) (\d+)$", log.read_text(), re.MULTILINE):
        counts[operation] += int(value)
    if set(counts) != operations or not all(value > 0 for value in counts.values()):
        raise ValueError(f"unexpected dispatch counts {dict(counts)}")
    if expected_calls is not None and any(value != expected_calls for value in counts.values()):
        raise ValueError(f"expected {expected_calls} calls to every boundary, got {dict(counts)}")
    print(json.dumps(dict(counts), indent=2))


if __name__ == "__main__":
    typer.run(main)
