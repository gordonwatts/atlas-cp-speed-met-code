"""Process and artifact receipts shared by the MET profiling commands.

Process-tree affinity observation follows the same /proc strategy as the previous
muon profiling harness. Linux wait4 measures the exact launched process tree.
"""
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import time

THREAD_ENV = {
    key:"1" for key in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS",
                        "NUMEXPR_NUM_THREADS","BLIS_NUM_THREADS","VECLIB_MAXIMUM_THREADS",
                        "ROOT_MAX_THREADS")
}
RELEASE_ENV = ("PYTHONPATH","PYTHONHOME","LD_LIBRARY_PATH","LD_PRELOAD","ROOTSYS",
               "ROOT_INCLUDE_PATH","CMAKE_PREFIX_PATH","ATLAS_LOCAL_ROOT_BASE","VIRTUAL_ENV",
               "AtlasVersion","AtlasProject","ROOT_VERSION","DATAPATH")


def json_write(path, data):
    with path.open("x", encoding="utf-8") as handle:
        json.dump(data, handle, indent=2, allow_nan=False)
        handle.write("\n")


def identify(path):
    path = path.resolve(strict=True)
    with path.open("rb") as handle:
        sha = hashlib.file_digest(handle, "sha256").hexdigest()
    return {"path":str(path), "sha256":sha, "bytes":path.stat().st_size}


def environment():
    env = os.environ.copy()
    for key in RELEASE_ENV:
        env.pop(key, None)
    env.update(THREAD_ENV)
    env.update(LC_ALL="C", PYTHONNOUSERSITE="1")
    return env


def process_tree(pid):
    seen, queue = set(), [pid]
    while queue:
        current = queue.pop()
        if current in seen:
            continue
        seen.add(current)
        proc = Path("/proc")/str(current)
        try:
            status = dict(line.split(":",1) for line in (proc/"status").read_text().splitlines() if ":" in line)
            if status["State"].strip().startswith("Z"):
                continue
            args = [os.fsdecode(item) for item in (proc/"cmdline").read_bytes().split(b"\0") if item]
            yield {"pid":current, "ppid":int(status["PPid"]), "affinity":status["Cpus_allowed_list"].strip(),
                   "threads":int(status["Threads"]), "argv":args, "exe":os.readlink(proc/"exe")}, proc
            for child_file in (proc/"task").glob("*/children"):
                try:
                    queue.extend(map(int, child_file.read_text().split()))
                except (FileNotFoundError, ProcessLookupError):
                    continue
        except (FileNotFoundError, ProcessLookupError):
            continue


def mapped_files(proc):
    result = set()
    try:
        for line in (proc/"maps").read_text().splitlines():
            fields = line.split(maxsplit=5)
            if len(fields) == 6 and fields[-1].startswith("/") and ".so" in Path(fields[-1]).name:
                result.add(fields[-1])
    except (FileNotFoundError, ProcessLookupError):
        pass
    return result


def launch(argv, out, demo, env):
    json_write(out/"command.json", argv)
    maps = set()
    seen = {}
    start = time.monotonic_ns()
    with (out/"stdout.log").open("x") as stdout, (out/"stderr.log").open("x") as stderr, (
            out/"processes.jsonl").open("x") as monitor:
        proc = subprocess.Popen(argv, stdout=stdout, stderr=stderr, env=env)
        while True:
            pid,status,usage = os.wait4(proc.pid,os.WNOHANG)
            if pid:
                proc.returncode = os.waitstatus_to_exitcode(status)
                break
            snapshot = []
            for row, root in process_tree(proc.pid):
                snapshot.append(row)
                if str(demo) in row["argv"] and Path(row["exe"]).name.startswith("python"):
                    identity = (row["affinity"],row["threads"])
                    seen.setdefault(row["pid"],set()).add(identity)
                    maps.update(mapped_files(root))
            monitor.write(json.dumps({"monotonic_ns":time.monotonic_ns(),"processes":snapshot})+"\n")
            monitor.flush()
            # Preserve elapsed finish without a fixed poll-interval bias. On
            # supported Linux this waits for the process FD up to the poll limit.
            if hasattr(os,"pidfd_open"):
                import select
                fd = os.pidfd_open(proc.pid)
                try:
                    with_poll = select.poll()
                    with_poll.register(fd, select.POLLIN)
                    with_poll.poll(200)
                finally:
                    os.close(fd)
            else:
                time.sleep(.2)
    end = time.monotonic_ns()
    receipt = {"wall_ns":end-start, "user_ns":round(usage.ru_utime*1e9),
               "system_ns":round(usage.ru_stime*1e9), "exit_code":proc.returncode,
               "started_monotonic_ns":start,"finished_monotonic_ns":end,
               "worker_observations":{str(key):[{"affinity":a,"threads":n} for a,n in sorted(val)]
                                      for key,val in sorted(seen.items())},
               "mapped_shared_objects":sorted(maps),
               "usage_scope":"Launched process plus waited descendants. Includes perf for sampled mode."}
    json_write(out/"process.json", receipt)
    if proc.returncode:
        raise RuntimeError(f"External workload exit {proc.returncode}. See {out/'stderr.log'}")
    return receipt


def counter_probe(perf, cpu, out, env):
    specs = {"software":("cpu-clock:u",),"hardware":("cycles:u","instructions:u")}
    evidence = {}
    for scope, events in specs.items():
        csv = out/f"{scope}.stat"
        command = ["taskset","-c",str(cpu),perf,"stat","-x",";","--no-big-num","-o",str(csv),
                   "-e",",".join(events),"--","/bin/true"]
        result = subprocess.run(command,text=True,capture_output=True,env=env)
        (out/f"{scope}.stderr").write_text(result.stderr)
        counts = {event:None for event in events}
        if csv.exists():
            for row in csv.read_text().splitlines():
                item = row.split(";")
                if len(item) < 3 or item[2].strip() not in counts:
                    continue
                try:
                    counts[item[2].strip()] = float(item[0]) if result.returncode == 0 else None
                except ValueError:
                    pass
        evidence[scope] = {"command":command,"exit":result.returncode,"numeric_returned_counters":counts}
    json_write(out/"counter_support.json",evidence)
    return evidence


def run_metadata(cpu, collector, python, demo, original, bundle, hash_input):
    manifest = json.loads(bundle.read_text())
    wheel = identify(bundle.parent/"wheels"/manifest["wheel"])
    if wheel["sha256"] != manifest["wheel_sha256"]:
        raise ValueError("Wheel hash differs from its bundle manifest")
    siblings = Path(f"/sys/devices/system/cpu/cpu{cpu}/topology/thread_siblings_list")
    stat = original.stat()
    footprint = {"path":str(original),"bytes":stat.st_size,"dev":stat.st_dev,"inode":stat.st_ino,
                 "mtime_ns":stat.st_mtime_ns}
    if hash_input:
        footprint["sha256"] = identify(original)["sha256"]
    return {"analysis_cpu":cpu,"controller_cpu":collector,"controller_effective_affinity":sorted(os.sched_getaffinity(0)),
            "analysis_cpu_thread_siblings":siblings.read_text().strip() if siblings.exists() else None,
            "platform":platform.platform(),"python":str(python),"demo":identify(demo),"wheel":wheel,
            "manifest":identify(bundle),"input":footprint,"thread_environment":THREAD_ENV,
            "removed_environment_names":RELEASE_ENV}
