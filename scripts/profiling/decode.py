"""Reduce perf-script observations with a fixed processing time and worker denominator.

Header and physical stack parsing are adapted from the reviewed small-R campaign
tools/native_stack_profile.py. This script reports MET-specific call anchors.
"""
from collections import Counter, defaultdict
from dataclasses import dataclass, field
import csv
from pathlib import Path
import re

HEADER = re.compile(r"^\s*(\d+)\s*/\s*(\d+)\s+\[(\d+)\]\s+(\d+)\.(\d+):\s+(\d+)\s+cpu-clock(?::u)?:\s*(.*)$")
FRAME = re.compile(r"^\s*(?:0x)?[\da-fA-F]+\s+(.*?)\s*\(([^()]+(?: \(deleted\))?)\)(?:\s+\(inlined\))?\s*$")
SUFFIX = re.compile(r"\+0x[\da-fA-F]+$")
UNKNOWN_NAMES = {"", "[unknown]", "unknown", "[unresolved]"}
UNKNOWN = ("[unknown]", "[unknown]")
LOSS = re.compile(r"\blost\s*[:=]?\s*(\d+)", re.IGNORECASE)

# Direct symbols identify a caller on a sample's own stack. Generic DSO, CPU
# or Python allocator ownership alone does not identify an analysis stage.
ANCHORS = (
    ("significance", ("atlas::met_significance::", "met_significance_event", "SignificanceNode")),
    ("track_soft_final", ("atlas::met::soft::", "SoftNode::", "FinalNode::")),
    ("maker", ("atlas::met::maker::", "MakerNode::", "maker::run", "MakerPlugin")),
    ("nnjvt", ("NnJvtUpdate::", "lwtnn::", "LWTNN::")),
    ("nominal_preparation", ("atlas::met::preparation::", "AfterCalibration::", "TrackRoles::",
                             "EgammaNode::", "TauNode::", "PointingProfile::", "LeakageGraph::")),
    ("shared_event_unattributed", ("met_columnar_full",)),
    ("histogram", ("boost::histogram", "pybind11_bh")),
    ("input_decompression", ("inflate_fast", "ZSTD_decompress", "LZ4_decompress", "lzma_code",
                             "cramjam::decompression")),
)


@dataclass
class Sample:
    pid: int
    tid: int
    cpu: int
    timestamp: int
    period: int
    stack: list = field(default_factory=list)


def frame(line):
    match = FRAME.fullmatch(line)
    if not match:
        return None
    symbol, dso = match.groups()
    name = SUFFIX.sub("", symbol).strip()
    return dso.strip(), ("[unknown]" if name in UNKNOWN_NAMES or name.startswith("0x") else name)


def records(stream, defects):
    active = None
    for line in stream:
        line = line.rstrip()
        header = HEADER.fullmatch(line)
        if header:
            if active:
                yield active
            pid, tid, cpu, sec, ns, period, inline = header.groups()
            if len(ns) > 9:
                raise ValueError("perf sample timestamp exceeds nanosecond precision")
            active = Sample(int(pid), int(tid), int(cpu), int(sec)*10**9 + int(ns.ljust(9,"0")), int(period))
            if inline:
                if first := frame(inline):
                    active.stack.append(first)
                else:
                    defects["unparsed_header_frames"] += 1
        elif "PERF_RECORD_LOST" in line:
            defects["lost_records"] += 1
            if count := LOSS.search(line):
                defects["lost_samples_reported"] += int(count[1])
            else:
                defects["lost_records_without_count"] += 1
        elif line.strip():
            if active and (entry := frame(line)):
                active.stack.append(entry)
            else:
                defects["unparsed_nonblank_lines"] += 1
    if active:
        yield active


def stage(stack):
    for dso, symbol in stack:
        for owner, patterns in ANCHORS:
            if any(anchor in symbol for anchor in patterns):
                return owner
        if "/boost_histogram/" in dso:
            return "histogram"
    return "stage_unattributed"


def rows(counter, denominator):
    for label, (samples, period) in sorted(counter.items(), key=lambda row:(-row[1][1], row[0])):
        yield (label, samples, period, 100. * period / denominator)


def table(path, header, entries):
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(entries)


def summarize(raw, out, worker, cpu, top):
    windows = [(item["start_monotonic_ns"], item["end_monotonic_ns"]) for item in worker["passes"]]
    if not windows or any(begin >= end for begin, end in windows):
        raise ValueError("Invalid recorded processing pass intervals")
    if any(a[1] > b[0] for a,b in zip(windows, windows[1:])):
        raise ValueError("Processing intervals overlap")

    def tally(group, key, period):
        group[key][0] += 1
        group[key][1] += period

    counters = Counter()
    defects = Counter()
    depths = Counter()
    self_cost, inclusive, paths, stage_cost = (defaultdict(lambda:[0,0]) for _ in range(4))
    total, count = 0, 0
    with raw.open(encoding="utf-8", errors="replace") as stream:
        for sample in records(stream, defects):
            in_window = any(begin <= sample.timestamp < end for begin,end in windows)
            counters["all_perf_samples"] += 1
            if sample.pid != worker["pid"]:
                counters["foreign_pid"] += 1
                continue
            if not in_window:
                counters["own_outside_processing"] += 1
                continue
            if sample.tid != worker["pid"] or sample.cpu != cpu:
                counters["unexpected_sample_tid_or_cpu_in_processing"] += 1
                continue
            if sample.period <= 0:
                counters["zero_period_in_processing"] += 1
                continue

            stack = sample.stack or [UNKNOWN]
            leaf = stack[0]
            period = sample.period
            count += 1
            total += period
            depths[sum(symbol != "[unknown]" for _,symbol in stack)] += 1
            tally(self_cost, leaf, period)
            tally(stage_cost, stage(stack), period)
            distinct = set(stack)
            for pair in distinct:
                tally(inclusive, pair, period)
            tally(paths, tuple(stack), period)

    if not count or not total:
        raise ValueError("No usable software-clock stack samples in the recorded processing windows")
    for group in (self_cost, stage_cost):
        if sum(value[0] for value in group.values()) != count or sum(value[1] for value in group.values()) != total:
            raise AssertionError("Exclusive cost tables do not reconcile with the whole profile")

    for cost, name in ((self_cost, "self"), (inclusive, "inclusive")):
        table(out / f"{name}.tsv", ("library","symbol","samples","period_ns","percent_of_processing"),
              ((*pair,samples,period,pct) for pair,samples,period,pct in rows(cost,total)))
    table(out / "stacks.tsv", ("stack_leaf_first","samples","period_ns","percent_of_processing"),
          ((" | ".join(f"{dso}::{symbol}" for dso,symbol in path),samples,period,pct)
           for path,samples,period,pct in rows(paths,total)))
    stage_rows = [{"stage":name,"samples":n,"period_ns":period,"percent_of_processing":pct}
                  for name,n,period,pct in rows(stage_cost,total)]
    table(out / "stage_by_stack.tsv", tuple(stage_rows[0]), (row.values() for row in stage_rows))

    resolved = list(row for row in rows(self_cost,total) if row[0][1] != "[unknown]")
    top_rows = [{"library":pair[0],"symbol":pair[1],"samples":n,"period_ns":period,"percent_of_processing":pct}
                for pair,n,period,pct in resolved[:top]]
    remainder = defaultdict(lambda:[0,0])
    for (dso,symbol),n,period,_ in resolved[top:]:
        remainder[dso][0] += n
        remainder[dso][1] += period
    remainder_rows = [{"library":dso,"samples":n,"period_ns":period,"percent_of_processing":pct}
                      for dso,n,period,pct in rows(remainder,total)]
    table(out/"named_remainder_by_library.tsv", ("library","samples","period_ns","percent_of_processing"),
          (row.values() for row in remainder_rows))
    unresolved_n = sum(n for (_,symbol),(n,_) in self_cost.items() if symbol == "[unknown]")
    unresolved_period = sum(period for (_,symbol),(_,period) in self_cost.items() if symbol == "[unknown]")
    if sum(row["period_ns"] for row in top_rows+remainder_rows) + unresolved_period != total:
        raise AssertionError("Top functions, named remainder and unresolved symbols do not sum")

    return {
        "denominator":{"event":"cpu-clock:u","sample_scope":"actual worker during recorded input passes",
                       "period_ns":total,"samples":count},
        "discarded_or_other":dict(counters),"decoder_diagnostics":dict(defects),
        "resolved_stack_depth_histogram":dict(sorted(depths.items())),
        "unresolved_leaf":{"samples":unresolved_n,"period_ns":unresolved_period,
                           "percent_of_processing":100.*unresolved_period/total},
        "top_self":top_rows,"stage_by_stack":stage_rows,"named_remainder_by_library":remainder_rows,
        "definitions":{"self":"Exclusive leaf share of recorded CPU-clock observations",
                       "inclusive":"Share of observations with the symbol anywhere in the stack, once per sample",
                       "stage":"First recognized semantic frame, or unattributed. Elapsed phase shares are independent."},
    }
