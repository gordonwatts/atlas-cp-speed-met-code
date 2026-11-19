#!/usr/bin/env python3
"""Measure the isolated faithful and M02 C endpoints in matched serial pairs."""
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import re
import statistics

import typer

# Run with the campaign scripts/ on PYTHONPATH. Share the faithful process and
# qualification rules instead of introducing a second timing instrument.
from benchmark_matched import NATIVE_C, PROBES, THREADS, quality, receipt
from benchmark_reference import run

SCOPES = ("full", "nominal")
DEFAULT_CHUNK_SIZE = 2000


def observe_variant(portable):
    bundle = portable/"bundle"
    wheels = list((bundle/"wheels").glob("*.whl"))
    if len(wheels)!=1:
        raise ValueError(f"Expected one wheel in {bundle}")
    python = portable/"pip-env"/"bin"/"python"
    demo = bundle/"met_demo.py"
    native = portable/"pip-env"/"lib"/"python3.14"/"site-packages"/"atlas_met_portable"/"native"
    identity = {str(path):receipt(path) for path in (
        python,demo,wheels[0],bundle/"bundle.json",portable/"installed-identities.json",
        portable/"pip-freeze.txt",*(native/name for name in NATIVE_C)
    )}
    expected = {str((native/name).resolve(strict=True)) for name in NATIVE_C}
    return python,demo,identity,expected


def main(input_file: Path, baseline: Path, trial: Path, stock_full: Path, stock_nominal: Path,
         output: Path, baseline_label: str = "faithful", trial_label: str = "m02",
         baseline_chunk_size: int | None = typer.Option(None,min=1),
         trial_chunk_size: int | None = typer.Option(None,min=1)):
    original = input_file.absolute()
    if not original.is_file():
        raise FileNotFoundError(original)
    labels = (baseline_label,trial_label)
    if len(set(labels)) != 2 or any(
        re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]*",label) is None or label=="mean_fractional_savings"
        for label in labels
    ):
        raise ValueError("Use distinct labels containing letters, digits, hyphens or underscores")
    variants = {name:observe_variant(path.resolve(strict=True))
                for name,path in zip(labels,(baseline,trial),strict=True)}
    chunk_size = dict(zip(labels,(baseline_chunk_size,trial_chunk_size),strict=True))
    effective = {role:requested or DEFAULT_CHUNK_SIZE for role,requested in chunk_size.items()}
    if variants[baseline_label][1]==variants[trial_label][1] and effective[baseline_label]==effective[trial_label]:
        raise ValueError("Use separate immutable bundles or distinct chunk sizes")
    expected_hist = {scope:receipt(path) for scope,path in zip(SCOPES,(stock_full,stock_nominal),strict=True)}
    standalone = {k:os.environ[k] for k in ("HOME","LANG","LC_ALL","TMPDIR") if k in os.environ}
    standalone.update(PATH="/usr/bin:/bin", **dict.fromkeys(THREADS,"1"))
    output = output.resolve()
    output.mkdir()
    proceedings = {}
    started = datetime.now(timezone.utc).isoformat()
    for round_index in range(4):
        scopes = SCOPES if round_index%2==0 else SCOPES[::-1]
        for scope in scopes:
            roles = tuple(variants) if round_index%2==0 else tuple(reversed(variants))
            for role in roles:
                dest = output/f"round-{round_index:02d}-{scope}-{role}"
                dest.mkdir()
                python,demo,_,expected = variants[role]
                histogram = dest/"histogram.npz"
                args = [str(python),str(demo),str(original),str(histogram),"--stop","20000"]
                if scope=="nominal":
                    args.append("--nominal-only")
                if chunk_size[role] is not None:
                    args.extend(("--chunk-size",str(chunk_size[role])))
                observed = run(args,standalone,dest,dest/"direct.log",PROBES)
                success = observed["exit_code"]==0 and histogram.exists()
                observed["native_identity_verified"] = set(observed["mapped_named_libraries"])==expected
                observed["histogram_verified"] = success and (
                    receipt(histogram)["sha256"]==expected_hist[scope]["sha256"])
                with (dest/"runtime.json").open("x") as handle:
                    json.dump(observed,handle,indent=2)
                    handle.write("\n")
                proceedings[scope,role,round_index] = observed

    summary = {}
    for scope in SCOPES:
        cohort = {}
        for role in variants:
            rows = [proceedings[scope,role,i] for i in range(4)]
            valid = all(row["valid"] and row["native_identity_verified"] and row["histogram_verified"]
                        for row in rows)
            cohort[role] = quality([row["wall_seconds"] for row in rows[1:]],valid)
            retained_rss = [row.get("peak_rss_kib") for row in rows[1:]]
            if all(isinstance(value,int) and value > 0 for value in retained_rss):
                cohort[role]["peak_rss_kib"] = {"retained":retained_rss,"maximum":max(retained_rss)}
        paired = [a-b for a,b in zip(cohort[baseline_label]["retained_seconds"],
                                     cohort[trial_label]["retained_seconds"],strict=True)]
        if all(row["qualified"] for row in cohort.values()):
            difference, spread = statistics.mean(paired), statistics.stdev(paired)
            cohort[f"paired_{baseline_label}_minus_{trial_label}_seconds"] = {
                "values":paired,"mean":difference,"sample_sd":spread}
            cohort["mean_fractional_savings"] = difference/cohort[baseline_label]["mean"]
        summary[scope] = cohort
    with (output/"summary.json").open("x") as handle:
        json.dump({"events":20000,"workload":"one direct-original histogram pass",
                   "reader_cpu":0,"excluded_round":0,"start_utc":started,
                   "variants":{"baseline":baseline_label,"trial":trial_label},
                   "requested_chunk_size":chunk_size,"effective_chunk_size":effective,
                   "end_utc":datetime.now(timezone.utc).isoformat(),
                   "identities":{role:variant[2] for role,variant in variants.items()},
                   "stock_histogram":expected_hist,"results":summary},handle,indent=2)
        handle.write("\n")
    if not all(summary[scope][role]["qualified"] for scope in SCOPES for role in variants):
        raise typer.Exit(1)


if __name__ == "__main__":
    typer.run(main)
