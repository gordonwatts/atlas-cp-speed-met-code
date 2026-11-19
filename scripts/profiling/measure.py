#!/usr/bin/env python3
"""Collect clean phase and native profiles of the accepted original-file MET C runner."""
from dataclasses import dataclass
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess

import numpy as np
import typer

from decode import summarize
from process import counter_probe, environment, identify, json_write, launch, run_metadata

app = typer.Typer(no_args_is_help=True)
UUID = "79ee6652-311f-11f1-aa37-0445010abeef"
ENTRIES = 20_000
STATES = ["NOSYS","MET_SoftTrk_ResoPara","MET_SoftTrk_ResoPerp","MET_SoftTrk_Scale__1down","MET_SoftTrk_Scale__1up"]
HIST_FIELDS = ["state_labels","bin_edges_gev","counts","sumw","sumw2"]
DIAG = ("import sys; sys.path.insert(0,sys.argv[1]); from phase_driver import run; "
        "run(sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5:])")
TOOLS = Path(__file__).resolve().parent


def read(path):
    return json.loads(path.read_text())


def app_summary(log):
    for line in reversed(log.read_text().splitlines()):
        try:
            value = json.loads(line)
        except json.JSONDecodeError:
            continue
        if isinstance(value, dict) and "histogram_sha256" in value:
            return value
    raise ValueError(f"Missing MET computation and pass summary in {log}")


def check_hist(candidate, expected, mode):
    labels = STATES if mode == "full" else STATES[:1]
    with np.load(candidate, allow_pickle=False) as observed, np.load(expected, allow_pickle=False) as stock:
        if observed.files != HIST_FIELDS or stock.files != HIST_FIELDS:
            raise ValueError("Expected the five accepted histogram fields in order")
        for field in HIST_FIELDS:
            a, b = observed[field], stock[field]
            if a.dtype != b.dtype or a.shape != b.shape or a.tobytes() != b.tobytes():
                raise ValueError(f"Fresh histogram differs from accepted stock field {field}")
        if observed["state_labels"].tolist() != labels or not np.all(observed["counts"].sum(axis=1) == ENTRIES):
            raise ValueError("Wrong accepted state or event count")
        digest = hashlib.sha256(observed["sumw"].tobytes()).hexdigest()
    return {"arrays_bitwise_equal":HIST_FIELDS,"saved_histogram":identify(candidate),
            "accepted_stock":identify(expected),"runner_sumw_sha256":digest}


def audit_loaded_paths(process, bundle, cpu):
    observations = process["worker_observations"]
    if len(observations) != 1 or next(iter(observations.values())) != [{"affinity":str(cpu),"threads":1}]:
        raise ValueError(f"Expected one observed CPython worker pinned to CPU {cpu}, got {observations}")
    manifest = read(bundle)
    expected = {Path(path).name: item["sha256"] for path,item in manifest["inputs"].items() if "/native/" in path}
    seen = {}
    for path in process["mapped_shared_objects"]:
        name = Path(path).name
        if "/cvmfs/" in path or name.startswith(("libCore.so","libRIO.so","libTree.so","libHist.so",
                                                "libMathCore.so","libxAOD","libAthena")):
            raise ValueError(f"Unpermitted release/ROOT library in installed workload: {path}")
        if name in expected:
            identity = identify(Path(path))
            if name in seen or identity["sha256"] != expected[name]:
                raise ValueError(f"Ambiguous or mismatched loaded runtime artifact {path}")
            seen[name] = identity
    if set(seen) != set(expected):
        raise ValueError(f"Process observer missed bundled native libraries {sorted(set(expected)-set(seen))}")
    return seen, int(next(iter(observations)))


def validate_phases(worker, process, repeats, digest, pid):
    if not worker["success"] or worker["pid"] != pid or len(worker["passes"]) != repeats:
        raise ValueError("Incomplete real-worker phase records")
    for index,row in enumerate(worker["passes"], 1):
        if row["pass"] != index or row["entries"] != ENTRIES or row["histogram_sha256"] != digest:
            raise ValueError("Repeat agreement or original pass visit record failed")
        if not 0 < row["cpu_ns"] or not 0 < row["wall_ns"]:
            raise ValueError("Invalid actual worker processing intervals")
    objects = worker["passes"][0]["objects"]
    if any(row["objects"] != objects for row in worker["passes"]):
        raise ValueError("Object visit counts changed in repeated original-file passes")
    script = worker["script_span"]
    phase_total = {unit:sum(x[unit] for x in worker["exclusive"].values()) for unit in script}
    residual = worker["unassigned_in_script"]
    external = {"wall_ns":process["wall_ns"],"cpu_ns":process["user_ns"]+process["system_ns"]}
    if any(script[unit] != phase_total[unit]+residual[unit] or residual[unit]<0 for unit in script):
        raise ValueError("Nonoverlapping phase budget does not reconcile")
    return {
        "original_objects":objects, "object_visits":{key:value*repeats for key,value in objects.items()},
        "processing_wall_ns":sum(row["wall_ns"] for row in worker["passes"]),
        "processing_cpu_ns":sum(row["cpu_ns"] for row in worker["passes"]),
        "worker_span":script, "worker_cpu_at_exec_finish_ns":worker["worker_cpu_at_exec_finish_ns"],
        "outside_worker_script_or_collector":{key:external[key]-script[key] for key in script},
        "unassigned_in_script":residual,
    }


@dataclass
class Workload:
    python: Path
    demo: Path
    original: Path
    bundle: Path
    cpu: int
    collector: int
    hash_input: bool
    chunk: int

    def initialize(self, out):
        if self.cpu == self.collector:
            raise ValueError("Choose a separate profiling controller logical CPU")
        os.sched_setaffinity(0,{self.collector})
        out.mkdir()
        info = run_metadata(self.cpu,self.collector,self.python,self.demo,self.original,self.bundle,self.hash_input)
        info["tools"] = {filename:identify(TOOLS/filename)
                         for filename in ("phase_driver.py","decode.py","measure.py","process.py")}
        info["chunk_size"] = self.chunk
        json_write(out/"provenance.json",info)

    def run(self, out, mode, repeats, expected, instrument=None, prefix=()):
        out.mkdir()
        hist, raw_phase = out/"histogram.npz",out/"phase_receipt.json"
        argv = [str(self.original),str(hist),"-n",str(repeats),"--chunk-size",str(self.chunk)]
        if mode == "nominal":
            argv.append("--nominal-only")
        command = ["taskset","-c",str(self.cpu),str(self.python)]
        command += (["-c",DIAG,str(TOOLS),str(self.demo),str(raw_phase),instrument]
                    if instrument else [str(self.demo)])
        process = launch([*prefix,*command,*argv],out,self.demo,environment())
        paths,pid = audit_loaded_paths(process,self.bundle,self.cpu)
        output = check_hist(hist,expected,mode)
        labels = STATES if mode == "full" else STATES[:1]
        app = app_summary(out/"stdout.log")
        if (app["entries"],app["event_visits"],app["repeats"],app["input_uuid"],app["states"],app["histogram_sha256"]) != (
                ENTRIES,repeats*ENTRIES,repeats,UUID,labels,output["runner_sumw_sha256"]):
            raise ValueError("Visible execution identity, original visits, or saved digest disagree")
        receipt = {"instrumentation":instrument,"runner_summary":app,"histogram":output,"loaded_libraries":paths,
                   "process":{key:process[key] for key in ("wall_ns","user_ns","system_ns","usage_scope")}}
        if instrument:
            worker = read(raw_phase)
            if worker["worker_affinity"] != [self.cpu] or worker["mode"] != instrument:
                raise ValueError("Unexpected worker instrumentation or CPU affinity")
            receipt["scope_accounting"] = validate_phases(worker,process,repeats,app["histogram_sha256"],pid)
            receipt["exclusive_phases"] = worker["exclusive"]
        json_write(out/"validated.json",receipt)
        return receipt


@app.callback()
def common(ctx: typer.Context, python: Path = typer.Option(...), demo: Path = typer.Option(...),
           input_file: Path = typer.Option(...,"--input"), bundle_manifest: Path = typer.Option(...),
           analysis_cpu: int = typer.Option(0,min=0), collector_cpu: int = typer.Option(1,min=0),
           chunk_size: int = typer.Option(2000,min=1), hash_input: bool = False):
    """Profile a frozen installed CLI in the same independent CPython environment."""
    interpreter = python.absolute()  # Resolving a venv interpreter symlink would bypass site-packages.
    if not interpreter.is_file() or not os.access(interpreter,os.X_OK):
        raise ValueError("Supply the tested portable venv interpreter")
    ctx.obj = Workload(interpreter,demo.resolve(strict=True),input_file.resolve(strict=True),
                       bundle_manifest.resolve(strict=True),analysis_cpu,collector_cpu,hash_input,chunk_size)


def ratio(slower, faster, key="wall_ns"):
    return slower["process"][key]/faster["process"][key]


@app.command()
def phases(ctx: typer.Context, out: Path = typer.Option(...), accepted_nominal: Path = typer.Option(...),
           accepted_full: Path = typer.Option(...)):
    """Measure named phase budgets for both single-pass modes plus direct control."""
    work = ctx.obj
    out = out.absolute()
    work.initialize(out)
    result = {}
    for mode,path in (("nominal",accepted_nominal),("full",accepted_full)):
        folder = out/mode
        folder.mkdir()
        control = work.run(folder/"control",mode,1,path.resolve(strict=True))
        detail = work.run(folder/"phase_timed",mode,1,path.resolve(strict=True),"phases")
        result[mode] = {"control":control,"phase_timed":detail,
                        "diagnostic_over_direct_wall_ratio":ratio(detail,control)}
    json_write(out/"summary.json",result)


@app.command()
def estimate(ctx: typer.Context, summary: Path, mode: str = typer.Option("full"),
             seconds: float = typer.Option(60.,min=60.), margin: float = typer.Option(1.20,min=1.)):
    """Choose repeats from the real phase-measured wall and process CPU time per pass."""
    binding = read(summary)[mode]["phase_timed"]["scope_accounting"]
    seconds_per_pass = min(binding["processing_wall_ns"],binding["processing_cpu_ns"])/1e9
    typer.echo(max(2,math.ceil(seconds*margin/seconds_per_pass)))


def perf_output(perf, out, filename, arguments):
    with (out/filename).open("x") as stdout, (out/(filename+".stderr")).open("x") as stderr:
        subprocess.run([perf,*arguments],stdout=stdout,stderr=stderr,env=environment(),check=True)


@app.command()
def native(ctx: typer.Context, out: Path = typer.Option(...), accepted: Path = typer.Option(...),
           repeats: int = typer.Option(...,"-n","--repeats",min=2), mode: str = typer.Option("full"),
           minimum_seconds: float = typer.Option(60.,min=60.), frequency: int = typer.Option(99,min=1,max=999),
           top: int = typer.Option(18,min=1), perf_path: str = "perf", pmu: bool = False):
    """Record native stacks across repeated true original-file passes and quantify harness overhead."""
    if mode not in ("nominal","full"):
        raise ValueError("Choose full or nominal")
    work = ctx.obj
    out, accepted = out.absolute(), accepted.resolve(strict=True)
    work.initialize(out)
    perf = shutil.which(perf_path)
    if not perf:
        json_write(out/"counter_support.json",{"perf":perf_path,"error":"Executable not found"})
        raise typer.Exit(2)
    json_write(out/"perf_tool.json",{"path":perf,
        "version":subprocess.check_output([perf,"--version"],text=True).strip(),
        "perf_event_paranoid":Path("/proc/sys/kernel/perf_event_paranoid").read_text().strip()})
    counters = counter_probe(perf,work.cpu,out,environment())
    if counters["software"]["numeric_returned_counters"]["cpu-clock:u"] is None:
        typer.echo(f"Software CPU counter unavailable. See {out/'counter_support.json'}",err=True)
        raise typer.Exit(2)

    direct = work.run(out/"control",mode,repeats,accepted)
    pass_timed = work.run(out/"pass_timed",mode,repeats,accepted,"passes")
    sample = out/"sampled"
    raw = sample/"perf.data"
    prefix = ["taskset","-c",str(work.collector),perf,"record","-e","cpu-clock:u","-F",str(frequency),
              "--call-graph","dwarf","--clockid","CLOCK_MONOTONIC","--timestamp","--period","--sample-cpu",
              "-o",str(raw),"--"]
    profiled = work.run(sample,mode,repeats,accepted,"passes",prefix)
    window = profiled["scope_accounting"]
    actual = min(window["processing_wall_ns"],window["processing_cpu_ns"])/1e9
    if actual < minimum_seconds:
        json_write(out/"insufficient_processing.json",{"observed_minimum_seconds":actual,"requested":minimum_seconds,
            "next_repeats":math.ceil(repeats*minimum_seconds/max(actual,.001)*1.15)})
        raise ValueError("Repeated original-file processing did not meet the coverage floor")
    raw_stacks = sample/"physical_stacks.txt"
    perf_output(perf,sample,raw_stacks.name,["script","--ns","--demangle","--no-inline","--show-lost-events",
                                           "-F","pid,tid,cpu,time,period,event,ip,sym,dso","-i",str(raw)])
    perf_output(perf,sample,"build_ids.txt",["buildid-list","-i",str(raw)])
    breakdown = summarize(raw_stacks,sample,read(sample/"phase_receipt.json"),work.cpu,top)
    if breakdown["discarded_or_other"].get("unexpected_sample_tid_or_cpu_in_processing",0):
        raise ValueError("Native observations during analysis used an unexpected TID or CPU")
    collector_rows = set()
    for line in (sample/"processes.jsonl").read_text().splitlines():
        for worker in json.loads(line)["processes"]:
            if Path(worker["exe"]).name == Path(perf).resolve().name:
                collector_rows.add((worker["pid"],worker["affinity"]))
    if not collector_rows or {affinity for pid,affinity in collector_rows} != {str(work.collector)}:
        raise ValueError("Unexpected or unobserved perf collector affinity")

    hardware = {"mode":"not_requested" if not pmu else "unsupported"}
    if pmu and all(value is not None for value in counters["hardware"]["numeric_returned_counters"].values()):
        pmu_dir = out/"hardware"
        args = ["taskset","-c",str(work.collector),perf,"stat","-x",";","--no-big-num","-o",
                str(pmu_dir/"counters.stat"),"-e","cycles:u,instructions:u","--"]
        hardware = {"mode":"collected","artifact":str(pmu_dir/"counters.stat"),
                    "workload":work.run(pmu_dir,mode,repeats,accepted,prefix=args)}
    summary = {"mode":mode,"sample_frequency_hz":frequency,"processing_floor_seconds":minimum_seconds,
               "profile_scope":"Repeated original-file windows only. Whole-process raw profile retained.",
               "control":direct,"pass_timed":pass_timed,"sampled":profiled,"native":breakdown,"hardware":hardware,
               "observed_perf_collector_pid_affinity":sorted(collector_rows),
               "overhead":{"pass_marker_over_direct_wall":ratio(pass_timed,direct),
                           "sampling_over_pass_marker_wall":ratio(profiled,pass_timed),
                           "sampled_over_pass_marker_worker_cpu":window["worker_cpu_at_exec_finish_ns"] /
                           pass_timed["scope_accounting"]["worker_cpu_at_exec_finish_ns"]}}
    json_write(out/"summary.json",summary)


if __name__ == "__main__":
    app()
