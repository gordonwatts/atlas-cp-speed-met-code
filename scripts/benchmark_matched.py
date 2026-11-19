"""Measure serial production A/B/C histogram endpoints and the independent read control."""
from __future__ import annotations

from datetime import datetime, timezone
from hashlib import sha256
import json
import os
from pathlib import Path
import statistics
import sys

import typer

from benchmark_reference import run

app = typer.Typer(add_completion=False)
NATIVE_C = (
    "libatlas_met_columnar.so", "libatlas_met_maker_correction.so",
    "libatlas_met_maker_provider.so", "libmet_preparation_provider.so", "libcorrectionlib.so",
)
FORBIDDEN_C = ("libCore.so", "libRIO.so", "libHist.so", "libTree.so", "libMathCore.so")
PROBES = (*NATIVE_C, *FORBIDDEN_C, "libAtlasMetFramework.so", "libMetAnalysisAlgorithmsLib.so")
SCOPES = ("nominal", "full", "stored")
THREADS = (
    "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS", "BLIS_NUM_THREADS", "RAYON_NUM_THREADS",
)


def receipt(path: Path) -> dict:
    resolved = path.resolve(strict=True)
    return {"path": str(resolved), "sha256": sha256(resolved.read_bytes()).hexdigest()}


def quality(values: list[float], process_ok: bool) -> dict:
    mean = statistics.mean(values)
    sd = statistics.stdev(values)
    median = statistics.median(values)
    variation = sd / mean
    largest = max(abs(value - median) / median for value in values)
    return {
        "retained_seconds": values, "mean": mean, "sample_sd": sd,
        "sd_over_mean": variation, "max_relative_median_distance": largest,
        "process_and_output_valid": process_ok,
        "qualified": process_ok and variation <= 0.05 and largest <= 0.10,
    }


def summarize(trials: list[dict]) -> dict:
    summary = {}
    for scope in SCOPES:
        chosen = [t for t in trials if t["scope"] == scope]
        if not chosen:
            continue
        scope_report = {}
        for path in dict.fromkeys(t["path"] for t in chosen):
            rows = sorted((t for t in chosen if t["path"] == path), key=lambda trial: trial["round"])
            scope_report[path] = {
                leg: quality(
                    [t[leg] for t in rows[1:]],
                    all(all(t[key] for key in observations) and t["histogram_verified"] for t in rows),
                )
                for leg, observations in (
                    ("framework_seconds", ("framework_valid",)),
                    ("histogram_seconds", ("histogram_valid",)),
                    ("total_seconds", ("framework_valid", "histogram_valid") if path != "C" else ("direct_valid",)),
                ) if leg in rows[0]
            }
        if "B" in scope_report:
            left, right = (scope_report[x]["framework_seconds"] for x in ("A", "B"))
            ratio = right["mean"] / left["mean"] if left["qualified"] and right["qualified"] else None
            scope_report["framework_ratio_B_over_A"] = ratio
            scope_report["framework_budget_1_20"] = ratio <= 1.20 if ratio is not None else None
        if "A" in scope_report and "C" in scope_report:
            left, right = (scope_report[x]["total_seconds"] for x in ("A", "C"))
            ratio = left["mean"] / right["mean"] if left["qualified"] and right["qualified"] else None
            scope_report["endpoint_ratio_A_over_C"] = ratio
        summary[scope] = scope_report
    return summary


@app.command()
def benchmark(
    input_file: Path, frozen_build: Path, portable: Path, output: Path, histogram_python: Path,
    label: str = "initial-matched", framework_budget_only: bool = False,
) -> None:
    repo = Path(__file__).resolve().parents[1]
    # Keep the user-facing .root alias. CPRun rejects its .art symlink target.
    input_file = input_file.absolute()
    if not input_file.is_file():
        raise FileNotFoundError(input_file)
    output = output.resolve()
    if os.environ.get("AtlasVersion") != "25.2.117" or os.environ.get("LD_PRELOAD"):
        raise RuntimeError("Start a clean stock AnalysisBase 25.2.117 release shell.")
    production = repo / "config" / "production.yaml"
    stored_config = repo / "config" / "stored-core.yaml"
    stock = Path(os.environ["AnalysisBase_DIR"]) / "bin" / "CPRun.py"
    boot = frozen_build / "source" / "components" / "framework" / "bridge_entry.py"
    bridge = frozen_build / "build" / "libAtlasMetFramework.so"
    provider = frozen_build / "build" / "maker_adapter" / "libatlas_met_maker_provider.so"
    core = repo / "results" / "core-exact-build-001" / "install" / "correctionlib" / "lib" / "libcorrectionlib.so"
    maker = repo / "results" / "maker-component-001" / "maker.json"
    soft = repo / "results" / "soft-component-exact-001" / "bundle" / "corrections.json"
    significance = repo / "results" / "significance-event-payload-001" / "met_significance_event.json"
    portable_python = portable / "pip-env" / "bin" / "python"
    demo = portable / "bundle" / "met_demo.py"
    native = portable / "pip-env" / "lib" / "python3.14" / "site-packages" / "atlas_met_portable" / "native"
    expected_native = {str((native / name).resolve(strict=True)) for name in NATIVE_C}
    expected_hashes = {scope: receipt(repo / "results" / f"reference-production-{scope}-002" / "histogram.npz")
                       for scope in SCOPES}
    files = [production, stored_config, stock, boot, bridge, provider, core, maker, soft, significance,
             demo, portable / "bundle" / "bundle.json", portable / "pip-freeze.txt",
             portable / "installed-identities.json", *sorted((portable / "bundle" / "wheels").glob("*.whl"))]
    identities = {path.name: receipt(path) for path in files}
    env_a = {key: value for key, value in os.environ.items() if not key.startswith("ATLAS_MET_")}
    env_a.update(dict.fromkeys(THREADS, "1"))
    env_b = {
        **env_a,
        "ATLAS_MET_FRAMEWORK_BRIDGE": str(bridge.resolve()),
        "ATLAS_MET_PROVIDER": str(provider.resolve()),
        "ATLAS_MET_MAKER_PAYLOAD": str(maker.resolve()),
        "ATLAS_MET_SOFT_PAYLOAD": str(soft.resolve()),
        "ATLAS_MET_SIGNIFICANCE_PAYLOAD": str(significance.resolve()),
        "ATLAS_MET_STAGE": "full",
    }
    # Standalone commands inherit a normal process environment plus the thread cap.
    standalone = {key: env_a[key] for key in ("LANG", "LC_ALL", "HOME", "TMPDIR", *THREADS) if key in env_a}
    standalone["PATH"] = "/usr/bin:/bin"
    output.mkdir()
    start = datetime.now(timezone.utc).isoformat()
    scopes = SCOPES[:2] if framework_budget_only else SCOPES
    trials = []
    for round_index in range(4):
        rotated_scopes = scopes[round_index % len(scopes):] + scopes[:round_index % len(scopes)]
        for scope in rotated_scopes:
            paths = ("A", "B", "C")
            if scope == "stored":
                paths = ("A", "C")
            if framework_budget_only:
                paths = ("A", "B")
            shift = round_index % len(paths)
            for path in paths[shift:] + paths[:shift]:
                cwd = output / f"round-{round_index:02d}-{scope}-{path}"
                cwd.mkdir()
                histogram_output = cwd / "histogram.npz"
                trial = {"round": round_index, "scope": scope, "path": path}
                if path == "C":
                    args = [str(portable_python), str(demo), str(input_file), str(histogram_output), "--stop", "20000"]
                    if scope == "nominal":
                        args.append("--nominal-only")
                    if scope == "stored":
                        args = [str(portable_python), str(repo / "scripts" / "columnar_stored_core.py"),
                                str(input_file), str(histogram_output)]
                    direct = run(args, standalone, cwd, cwd / "direct.log", PROBES)
                    maps = set(direct["mapped_named_libraries"])
                    expected = expected_native if scope != "stored" else set()
                    exact = maps == expected
                    trial.update(direct=direct, direct_valid=direct["valid"] and exact,
                                 total_seconds=direct["wall_seconds"])
                    identity_ok = exact
                    success = direct["exit_code"] == 0
                else:
                    environment = env_a if path == "A" else env_b
                    runner = stock if path == "A" else boot
                    config = stored_config if scope == "stored" else production
                    args = [sys.executable, str(runner), "--direct-driver", "-i", str(input_file),
                            "-t", str(config), "-o", "reference", "-e", "20000"]
                    if scope != "full":
                        args.append("--no-systematics")
                    framework = run(args, environment, cwd, cwd / "framework.log", PROBES)
                    mappings = set(framework["mapped_named_libraries"])
                    expected = {str(x.resolve()) for x in (bridge, provider, core)}
                    verified = expected <= mappings if path == "B" else not expected & mappings
                    trial.update(framework=framework, framework_valid=framework["valid"] and verified,
                                 framework_seconds=framework["wall_seconds"])
                    identity_ok = verified
                    success = framework["exit_code"] == 0
                    if success:
                        hist_args = [str(histogram_python), str(repo / "scripts" / "histogram_reference.py"),
                                     "reference.root", "hist-reference.root", str(histogram_output)]
                        if scope == "stored":
                            hist_args.append("--stored")
                        histogram = run(hist_args, standalone, cwd, cwd / "histogram.log")
                        trial.update(histogram=histogram, histogram_valid=histogram["valid"],
                                     histogram_seconds=histogram["wall_seconds"],
                                     total_seconds=framework["wall_seconds"] + histogram["wall_seconds"])
                        success = histogram["exit_code"] == 0
                artifact = receipt(histogram_output) if histogram_output.is_file() else None
                trial.update(identity_verified=identity_ok, histogram_output=artifact,
                             histogram_verified=artifact is not None and artifact["sha256"] == expected_hashes[scope]["sha256"])
                with (cwd / "runtime.json").open("x") as log:
                    json.dump(trial, log, indent=2)
                    log.write("\n")
                if not success or not trial["histogram_verified"]:
                    raise RuntimeError(f"Incomplete or numerically changed timed endpoint {cwd}.")
                trials.append(trial)
                print(round_index, scope, path, f'{trial["total_seconds"]:.3f} sec', flush=True)
    report = {"label": label, "start_utc": start, "end_utc": datetime.now(timezone.utc).isoformat(),
              "events": 20000, "cpu": 0, "first_round_excluded": True,
              "identities": identities, "histogram_reference": expected_hashes, "trials": trials,
              "cohorts": summarize(trials)}
    with (output / "benchmark.json").open("x") as log:
        json.dump(report, log, indent=2)
        log.write("\n")
    with (output / "summary.json").open("x") as log:
        json.dump({key: value for key, value in report.items() if key != "trials"}, log, indent=2)
        log.write("\n")


if __name__ == "__main__":
    app()
