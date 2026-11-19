"""Measure sequential, matched stock A framework and histogram process endpoints."""
from __future__ import annotations

from datetime import datetime, timezone
import json
import os
from pathlib import Path
import resource
import statistics
import subprocess
import sys
import threading
import time

import typer

app = typer.Typer(add_completion=False)
CPU = 0
SCOPES = ("stored", "nominal", "full")


def snapshot_cpu() -> list[int]:
    line = next(line for line in Path("/proc/stat").read_text().splitlines() if line.startswith(f"cpu{CPU} "))
    return [int(value) for value in line.split()[1:]]


def cpu_seconds() -> float:
    usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    return usage.ru_stime + usage.ru_utime


def monitor(pid: int, done: threading.Event, rows: list) -> None:
    other_cpus = set(os.sched_getaffinity(0)) - {CPU}
    if other_cpus:
        os.sched_setaffinity(0, {min(other_cpus)})
    proc = Path(f"/proc/{pid}")
    while not done.is_set():
        try:
            status = dict(line.split(":", 1) for line in (proc / "status").read_text().splitlines() if ":" in line)
            executable = str((proc / "exe").readlink())
            children = (proc / "task" / str(pid) / "children").read_text().split()
            if "python" in Path(executable).name:
                rows.append({
                    "affinity": status["Cpus_allowed_list"].strip(),
                    "threads": int(status["Threads"]),
                    "children": children,
                    "executable": executable,
                })
        except (FileNotFoundError, ProcessLookupError, KeyError):
            break
        done.wait(0.1)


def run(command: list[str], environment: dict, cwd: Path, log: Path) -> dict:
    command = ["taskset", "-c", str(CPU), *command]
    events = []
    finished = threading.Event()
    with log.open("x") as output:
        before_host = snapshot_cpu()
        before_cpu = cpu_seconds()
        start = time.perf_counter_ns()
        child = subprocess.Popen(command, cwd=cwd, env=environment, stdout=output, stderr=subprocess.STDOUT)
        sampler = threading.Thread(target=monitor, args=(child.pid, finished, events))
        sampler.start()
        status = child.wait()
        end = time.perf_counter_ns()
        finished.set()
        sampler.join()
        process_cpu = cpu_seconds() - before_cpu
        after_host = snapshot_cpu()
    observation = {
        "arguments": command,
        "exit_code": status,
        "wall_seconds": (end - start) * 1e-9,
        "process_cpu_seconds": process_cpu,
        "host_cpu0_delta_ticks": [after - before for before, after in zip(before_host, after_host, strict=True)],
        "samples": events,
    }
    if status != 0:
        observation["valid"] = False
    else:
        observation["valid"] = bool(events) and all(
            row["affinity"] == str(CPU) and row["threads"] == 1 and not row["children"] for row in events
        )
    return observation


def cohort(results: list[dict]) -> dict:
    metrics = {}
    for scope in SCOPES:
        selected = [trial for trial in results if trial["scope"] == scope]
        metrics[scope] = {}
        for leg in ("framework_seconds", "histogram_seconds", "total_seconds"):
            values = [trial[leg] for trial in selected[1:]]
            mean = statistics.mean(values)
            sd = statistics.stdev(values)
            median = statistics.median(values)
            maximum_deviation = max(abs(value - median) / median for value in values)
            metrics[scope][leg] = {
                "retained_seconds": values,
                "mean": mean,
                "sample_sd": sd,
                "sd_over_mean": sd / mean,
                "maximum_fractional_deviation_from_median": maximum_deviation,
                "qualified": sd / mean <= 0.05 and maximum_deviation <= 0.1
                and all(trial["framework"]["valid"] and trial["histogram"]["valid"] for trial in selected),
            }
    return metrics


@app.command()
def benchmark(
    input_file: Path,
    output: Path,
    histogram_python: Path,
    label: str = "initial",
    events: int = typer.Option(20000, min=1),
) -> None:
    repo = Path(__file__).resolve().parents[1]
    input_file = input_file.absolute()
    output = output.resolve()
    histogram_python = histogram_python.absolute()
    output.mkdir()
    if os.environ.get("AtlasVersion") != "25.2.117":
        raise ValueError("Source the stock AnalysisBase 25.2.117 environment before benchmarking.")
    if os.environ.get("LD_PRELOAD"):
        raise ValueError("Benchmark with a clean stock release process environment.")
    release_python = sys.executable
    cprun = str(Path(os.environ["AnalysisBase_DIR"]) / "bin" / "CPRun.py")
    framework_env = dict(os.environ)
    histogram_env = dict(os.environ)
    for key in ("PYTHONHOME", "PYTHONPATH", "LD_LIBRARY_PATH"):
        histogram_env.pop(key, None)
    start = datetime.now(timezone.utc).isoformat()
    results = []
    for round_index in range(4):
        ordered_scopes = SCOPES[round_index % len(SCOPES):] + SCOPES[:round_index % len(SCOPES)]
        for scope in ordered_scopes:
            cwd = output / f"round-{round_index:02d}-{scope}"
            cwd.mkdir()
            config = repo / "config" / ("stored-core.yaml" if scope == "stored" else "production.yaml")
            args = [release_python, cprun, "--direct-driver", "-i", str(input_file), "-t", str(config), "-o", "reference", "-e", str(events)]
            if scope != "full":
                args.append("--no-systematics")
            framework = run(args, framework_env, cwd, cwd / "framework.log")
            trial = {"scope": scope, "round": round_index, "framework": framework}
            if framework["exit_code"] == 0:
                hist_args = [
                    str(histogram_python), str(repo / "scripts" / "histogram_reference.py"),
                    "reference.root", "hist-reference.root", "histogram.npz",
                ]
                if scope == "stored":
                    hist_args.append("--stored")
                histogram = run(hist_args, histogram_env, cwd, cwd / "histogram.log")
                trial["histogram"] = histogram
                trial["framework_seconds"] = framework["wall_seconds"]
                trial["histogram_seconds"] = histogram["wall_seconds"]
                trial["total_seconds"] = framework["wall_seconds"] + histogram["wall_seconds"]
            with (cwd / "runtime.json").open("x") as report:
                json.dump(trial, report, indent=2)
                report.write("\n")
            if framework["exit_code"] != 0 or trial["histogram"]["exit_code"] != 0:
                raise RuntimeError(f"Endpoint failed in {cwd}.")
            results.append(trial)
            print(round_index, scope, f"{trial['total_seconds']:.3f} sec", flush=True)
    with (output / "benchmark.json").open("x") as file:
        json.dump(
            {
                "label": label, "start_utc": start, "events": events,
                "meters": "Fresh subprocess starts to exit, separate framework and histogram legs.",
                "trials": results, "cohorts": cohort(results),
            },
            file, indent=2,
        )
        file.write("\n")


if __name__ == "__main__":
    app()
