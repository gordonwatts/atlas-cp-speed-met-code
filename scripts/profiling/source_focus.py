"""Observe the original preparation projection with an exact saved-histogram control."""
from __future__ import annotations

import ast
from collections import Counter, defaultdict
from contextlib import ExitStack, contextmanager
import cProfile
from functools import wraps
import inspect
import json
from pathlib import Path
import pstats
import runpy
import sys
from textwrap import dedent
import time
from unittest.mock import patch

import atlas_met_portable  # Load the wheel's runtime before native Python dependencies.
import awkward as ak
import typer
import uproot

import atlas_met_columnar as transport
import atlas_met_prepare as prep

app = typer.Typer(add_completion=False)


class Probe:
    def __init__(self):
        self.scope = []
        self.counts = defaultdict(lambda: defaultdict(Counter))
        self.locations = Counter()
        self.raw_array_uses = Counter()
        self.repeated_cluster_accesses = Counter()

    @contextmanager
    def observe(self, group, name, rows=0):
        table = self.counts[group][name]
        table["calls"] += 1
        table["rows"] += rows
        start = time.perf_counter_ns()
        try:
            yield
        finally:
            table["elapsed_ns"] += time.perf_counter_ns() - start

    @contextmanager
    def species(self, name):
        self.scope.append(name)
        try:
            with self.observe("species", name):
                yield
        finally:
            self.scope.pop()

    def current(self):
        return self.scope[-1] if self.scope else "outside_object_preparation"

    @contextmanager
    def candidate(self, projection, method, field, target, rows):
        name = target + "AuxDyn." + field
        available = name in projection.data.raw
        kind = "absent_field" if not available else "nonempty_rows" if len(rows) else "empty_rows"
        scope = self.current()
        if method == "sampling" and available:
            self.repeated_cluster_accesses[(scope, str(name))] += 1
        if available:
            self.raw_array_uses[(scope, method, name, kind)] += 1
        with self.observe("candidate_group", f"{scope}/{method}/{kind}", len(rows)):
            with self.observe("candidate_detail", f"{scope}/{method}/{kind}/{name}", len(rows)):
                yield


def targeted_method(method, probe):
    source = dedent(inspect.getsource(getattr(prep.Resolved, method)))
    tree = ast.parse(source)
    function, = tree.body
    loop, = (row for row in function.body if isinstance(row, ast.For))
    if ast.unparse(loop.iter) != "self.locations" or ast.unparse(loop.target) != "(target, local, rows)":
        raise ValueError(f"Unexpected accepted candidate loop in Resolved.{method}")
    marker = ast.parse(f'with _candidate_probe(self, "{method}", field, target, rows): pass').body[0]
    marker.body = loop.body
    loop.body = [marker]
    ast.fix_missing_locations(tree)
    namespace = {**vars(prep), "_candidate_probe": probe.candidate}
    exec(compile(tree, f"<focused Resolved.{method}>", "exec"), namespace)
    return namespace[method]


def tracked(probe, group, label, function):
    @wraps(function)
    def capture(*args, **kwargs):
        with probe.observe(group, label):
            return function(*args, **kwargs)
    return capture


def run_audit(demo, input_file, output, probe):
    def species(function, name=None):
        @wraps(function)
        def capture(*args, **kwargs):
            with probe.species(name or args[2]):
                return function(*args, **kwargs)
        return capture

    original_init = prep.Resolved.__init__

    def resolved_init(instance, *args, **kwargs):
        original_init(instance, *args, **kwargs)
        for target, _, rows in instance.locations:
            probe.locations[(probe.current(), target, "empty" if len(rows) == 0 else "nonempty")] += 1

    arrays = uproot.behaviors.TBranch.HasBranches.arrays

    @wraps(arrays)
    def observed_arrays(owner, *args, **kwargs):
        role = "plain_collection_branches" if owner.name == "CollectionTree" else "nested_owner_branches"
        with probe.observe("read", role):
            return arrays(owner, *args, **kwargs)

    operations = ("flatten", "pad_none", "is_none", "fill_none", "to_numpy", "num", "to_packed")
    with ExitStack() as stack:
        for name in ("column", "link_fields", "sampling"):
            stack.enter_context(patch.object(prep.Resolved, name, targeted_method(name, probe)))
        stack.enter_context(patch.object(prep.Resolved, "__init__", resolved_init))
        stack.enter_context(patch.object(prep, "prepare_muon", species(prep.prepare_muon, "muon")))
        stack.enter_context(patch.object(prep, "prepare_tau", species(prep.prepare_tau, "tau")))
        stack.enter_context(patch.object(prep, "prepare_egamma", species(prep.prepare_egamma)))
        for operation in operations:
            original = getattr(ak, operation)

            def capture(*args, name=operation, original=original, **kwargs):
                caller = sys._getframe(1).f_code.co_name
                with probe.observe("awkward", f"{probe.current()}/{caller}/{name}"):
                    return original(*args, **kwargs)

            stack.enter_context(patch.object(ak, operation, capture))
        for method in ("resolve", "own", "offset"):
            stack.enter_context(patch.object(
                prep.Data, method, tracked(probe, "link_and_gather", method, getattr(prep.Data, method))
            ))
        stack.enter_context(patch.object(uproot.behaviors.TBranch.HasBranches, "arrays", observed_arrays))
        stack.enter_context(patch.object(
            transport, "read_associations",
            tracked(probe, "read", "nested_association_projection", transport.read_associations)
        ))
        stack.enter_context(patch.object(
            transport.MemberwiseJets, "read",
            tracked(probe, "read", "memberwise_jet_link_projection", transport.MemberwiseJets.read)
        ))
        prior = sys.argv
        try:
            sys.argv = [str(demo), str(input_file), str(output), "--stop", "20000"]
            try:
                runpy.run_path(str(demo), run_name="__main__")
            except SystemExit as exit:
                if exit.code:
                    raise
        finally:
            sys.argv = prior


@app.command()
def main(demo: Path, input_file: Path, accepted_histogram: Path, out: Path):
    out.mkdir()
    probe = Probe()
    profile = cProfile.Profile()
    output = out/"histogram.npz"
    with profile:
        run_audit(demo, input_file, output, probe)
    if output.read_bytes() != accepted_histogram.read_bytes():
        raise ValueError("Diagnostic output disagrees with the accepted stock histogram")
    profile.dump_stats(out/"python-calls.prof")
    with (out/"python-calls.txt").open("x") as handle:
        stats = pstats.Stats(profile, stream=handle).strip_dirs().sort_stats("cumulative")
        stats.print_stats(45)
        stats.print_callers("flatten|to_packed|pad_none|sampling|resolve|<focused|basket|decode", 65)
    receipt = {
        "events": 20000,
        "histogram_identical_to_stock": True,
        "meters": "Instrumented scopes overlap. Use their times for direction and the application phase control for whole-pass elapsed.",
        "resolved_init_locations": {"/".join(key): value for key, value in sorted(probe.locations.items())},
        "resolved_field_accesses": {"/".join(key): value for key, value in sorted(probe.raw_array_uses.items())},
        "calorimeter_sampling_accesses": {"/".join(key): value for key, value in sorted(probe.repeated_cluster_accesses.items())},
        "diagnostic_scopes": probe.counts,
    }
    with (out/"source-focus.json").open("x") as file:
        json.dump(receipt, file, indent=2)
        file.write("\n")


if __name__ == "__main__":
    app()
