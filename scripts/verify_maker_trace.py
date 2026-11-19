"""Compare ordered B maker decisions with the independent stock trace."""
from __future__ import annotations

from collections import Counter
from itertools import zip_longest
import json
from pathlib import Path

import typer

app = typer.Typer(add_completion=False)
CHECKED_ACTIONS = {"associationFlag", "jetJvt"}


def checked_rows(path: Path, entries: int):
    with path.open() as source:
        for line in source:
            entry, state, action, *data = line.rstrip("\n").split("\t")
            if int(entry) >= entries:
                break
            if action in CHECKED_ACTIONS:
                yield [int(entry), int(state), action, *map(int, data)]


@app.command()
def compare(stock: Path, candidate: Path, entries: int, output: Path) -> None:
    counts = Counter()
    first = None
    for a, b in zip_longest(checked_rows(stock, entries), checked_rows(candidate, entries)):
        if a != b:
            first = {"stock": a, "candidate": b}
            break
        counts[a[2]] += 1
    passed = first is None and all(counts[category] > 0 for category in CHECKED_ACTIONS)
    with output.open("x") as file:
        json.dump({"entries": entries, "exact_ordered_counts": counts, "first_difference": first, "passed": passed}, file, indent=2)
        file.write("\n")
    if not passed:
        raise typer.Exit(1)


if __name__ == "__main__":
    app()
