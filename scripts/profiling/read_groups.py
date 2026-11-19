"""Measure TBasket materialization and value interpretation by original ROOT branch."""
from __future__ import annotations

from collections import Counter, defaultdict
from contextlib import ExitStack
from functools import wraps
from hashlib import sha256
import json
from pathlib import Path
import runpy
import sys
import time
from unittest.mock import patch

import atlas_met_portable  # Establish the bundled runtime before uproot.
import typer
import uproot

import atlas_met_prepare as preparation

app = typer.Typer(add_completion=False)
INTERPRETATIONS = (
    uproot.interpretation.numerical.AsDtype,
    uproot.interpretation.jagged.AsJagged,
    uproot.interpretation.objects.AsObjects,
    uproot.interpretation.strings.AsStrings,
)
TBASKET = uproot.models.TBasket.Model_TBasket


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


@app.command()
def measure(demo: Path, input_file: Path, accepted_histogram: Path, output: Path):
    output.mkdir()
    histogram = output / "histogram.npz"
    leaf = defaultdict(Counter)
    process = Counter()
    read_members = TBASKET.read_members
    prepare_chunk = preparation.prepare_chunk
    depth = 0

    @wraps(read_members)
    def read(basket, chunk, cursor, context, file):
        start = time.process_time_ns()
        result = read_members(basket, chunk, cursor, context, file)
        spent = time.process_time_ns() - start
        if context.get("read_basket", True):
            stats = leaf[basket._parent.name]
            stats["basket_materializations"] += 1
            stats["basket_materialization_cpu_ns"] += spent
            stats["compressed_bytes_visited_with_key"] += basket.compressed_bytes
            stats["uncompressed_bytes_visited_with_key"] += basket.uncompressed_bytes
        return result

    def interpreted(function):
        @wraps(function)
        def decode(interpretation, data, byte_offsets, basket, branch, *args, **kwargs):
            nonlocal depth
            depth += 1
            if depth != 1:
                try:
                    return function(interpretation, data, byte_offsets, basket, branch, *args, **kwargs)
                finally:
                    depth -= 1
            start = time.process_time_ns()
            try:
                return function(interpretation, data, byte_offsets, basket, branch, *args, **kwargs)
            finally:
                spent = time.process_time_ns() - start
                depth -= 1
                stats = leaf[branch.name]
                stats["outer_interpretations"] += 1
                stats["outer_interpretation_cpu_ns"] += spent
        return decode

    @wraps(prepare_chunk)
    def prepare(*args, **kwargs):
        start = time.process_time_ns()
        try:
            return prepare_chunk(*args, **kwargs)
        finally:
            process["nominal_preparation_calls"] += 1
            process["nominal_preparation_cpu_ns"] += time.process_time_ns() - start

    prior = sys.argv
    wall_start = time.perf_counter_ns()
    cpu_start = time.process_time_ns()
    with ExitStack() as stack:
        stack.enter_context(patch.object(TBASKET, "read_members", read))
        stack.enter_context(patch.object(preparation, "prepare_chunk", prepare))
        for interpretation in INTERPRETATIONS:
            stack.enter_context(patch.object(interpretation, "basket_array", interpreted(interpretation.basket_array)))
        try:
            sys.argv = [str(demo), str(input_file), str(histogram), "--stop", "20000"]
            try:
                runpy.run_path(str(demo), run_name="__main__")
            except SystemExit as exit:
                if exit.code:
                    raise
        finally:
            sys.argv = prior
    process["diagnostic_launcher_and_execution_cpu_ns"] = time.process_time_ns() - cpu_start
    process["diagnostic_launcher_and_execution_wall_ns"] = time.perf_counter_ns() - wall_start
    if histogram.read_bytes() != accepted_histogram.read_bytes() or not leaf or not process["nominal_preparation_calls"]:
        raise ValueError("No exact full original-file evidence")

    groups = defaultdict(Counter)
    for branch, stats in leaf.items():
        head = branch.split("AuxDyn.", 1)[0]
        groups[head].update(stats)
    with (output / "read-groups.json").open("x") as result:
        json.dump({
            "events": 20000,
            "script_sha256": digest(demo),
            "histogram_sha256": digest(histogram),
            "notes": "Process CPU in named diagnostic hooks. Basket reads include ROOT key parsing and decompression; interpretation counts outer value conversion. Visits include duplicates. Whole-phase residual includes branch dispatch, packing and array finalization. Use logical ROOT footprint for distinct stored byte extents.",
            "process": process,
            "branch_groups": dict(sorted(groups.items())),
            "branches": dict(sorted(leaf.items())),
        }, result, indent=2)
        result.write("\n")


if __name__ == "__main__":
    app()
