#!/usr/bin/env python3
"""Render compact technical tables from MET diagnostic receipts."""
import json
from pathlib import Path
from typing import Annotated

import typer

PHASE_LABELS = {
    "imports":"Python and package imports",
    "initialize":"Resolve payloads and initialize native context",
    "pass_reset":"Reset per-pass native state",
    "book_histogram":"Book histogram",
    "open_original":"Open original PHYSLITE",
    "input_metadata":"Select branches and initialize association reader",
    "input_read_decode":"Read and decode original branches and associations",
    "raw_transport":"Pack raw typed buffers",
    "nominal_projection":"Project original nominal inputs and links",
    "shared_event_evaluation":"Evaluate shared MET event and state nodes",
    "fill_histogram":"Fill histogram",
    "pass_histogram_signature":"Compute repeat agreement digest",
    "close_native":"Close native context",
    "save_histogram":"Validate and save histogram",
}


def ms(value):
    return f"{value/1e6:,.2f}"


def sec(value):
    return f"{value/1e9:,.3f}"


def line(value):
    print(value)


def table(header, rows):
    line("| " + " | ".join(header) + " |")
    line("| " + " | ".join("---" for _ in header) + " |")
    for row in rows:
        cells = [str(value).replace("|", "\\|").replace("\n", " ") for value in row]
        line("| " + " | ".join(cells) + " |")
    line("")


def phase_table(mode, data):
    instrumented = data["phase_timed"]
    scope = instrumented["scope_accounting"]
    control = data["control"]["process"]
    line(f"## {mode.capitalize()} phase accounting")
    line("")
    line("Elapsed and process CPU values are for one complete original-file pass. "
         "Phase rows are exclusive within the diagnostic worker. The control runs the visible script directly.")
    line("")
    lines = []
    native = [name for name in instrumented["exclusive_phases"] if name.startswith("nominal_node.")]
    for key,label in PHASE_LABELS.items():
        row = instrumented["exclusive_phases"].get(key)
        if row:
            lines.append((label, ms(row["wall_ns"]), ms(row["cpu_ns"]),row["calls"]))
        if key == "nominal_projection":
            for name in native:
                sub = instrumented["exclusive_phases"][name]
                lines.append((f"Native nominal call {name.removeprefix('nominal_node.')}",
                              ms(sub["wall_ns"]),ms(sub["cpu_ns"]),sub["calls"]))
    omitted = set(instrumented["exclusive_phases"]) - set(PHASE_LABELS) - set(native)
    for key in sorted(omitted):
        row = instrumented["exclusive_phases"][key]
        lines.append((key,ms(row["wall_ns"]),ms(row["cpu_ns"]),row["calls"]))
    for key,label in (("unassigned_in_script","Unassigned inside visible script"),
                      ("worker_span","Diagnostic visible script total"),
                      ("outside_worker_script_or_collector","Diagnostic subprocess outside visible script")):
        row = scope[key]
        lines.append((label,ms(row["wall_ns"]),ms(row["cpu_ns"]),""))
    lines.append(("Direct subprocess total",ms(control["wall_ns"]),ms(control["user_ns"]+control["system_ns"]),""))
    table(("Phase","Elapsed ms","Process CPU ms","Calls"),lines)


def native_tables(summary):
    shard = summary["native"]
    scope = summary["sampled"]["scope_accounting"]
    app = summary["sampled"]["runner_summary"]
    den = shard["denominator"]
    line("## Native stack observations")
    line("")
    line(f"The {summary['mode']} worker completed {app['repeats']} original-file passes "
         f"over {app['entries']:,} unique events and {app['event_visits']:,} event visits. "
         f"Recorded processing windows total {sec(scope['processing_wall_ns'])} s elapsed "
         f"and {sec(scope['processing_cpu_ns'])} s process CPU.")
    line("")
    table(("Object","Original rows","Profile visits"),
          ((k, f"{v:,}", f"{scope['object_visits'][k]:,}") for k,v in scope["original_objects"].items()))
    line(f"The native denominator is {den['samples']:,} processing-window CPU-clock samples "
         f"weighted by {den['period_ns']:,} ns of recorded periods. "
         "Exclusive share uses the executing leaf. Inclusive tables count a symbol wherever it appears on a stack "
         "and overlap across symbols. Sample shares describe observed CPU activity.")
    line("")
    rows = [(r["library"],r["symbol"],r["samples"],f"{r['percent_of_processing']:.2f}%")
            for r in shard["top_self"]]
    missing = shard["unresolved_leaf"]
    rows.append(("(any)","Unresolved leaf",missing["samples"],f"{missing['percent_of_processing']:.2f}%"))
    table(("Library","Top self symbol","Samples","Exclusive CPU share"),rows)
    table(("Library","Remaining named leaf samples","Exclusive CPU share"),
          ((r["library"],r["samples"],f"{r['percent_of_processing']:.2f}%") for r in shard["named_remainder_by_library"]))
    table(("Recognized stack operation","Samples","Exclusive allocated sample share"),
          ((r["stage"],r["samples"],f"{r['percent_of_processing']:.2f}%") for r in shard["stage_by_stack"]))
    line("Stage assignment uses explicit native stack anchors. The unattributed stage includes known leaf functions "
         "whose workflow caller could not be identified. The unresolved-leaf row is a separate symbol-coverage measure.")
    line("")
    line(f"Perf diagnostics {json.dumps(shard['decoder_diagnostics'], sort_keys=True)}. "
         f"Resolved frame depths {json.dumps(shard['resolved_stack_depth_histogram'], sort_keys=True)}. "
         f"Hardware counters {summary['hardware']['mode']}.")
    line("")
    table(("Matched diagnostic overhead ratio","Ratio"),
          ((key.replace("_"," "),f"{value:.3f}") for key,value in summary["overhead"].items()))
    line("The profiler observes the sustained repeated-input workload. The direct single-pass benchmark has its own timing ledger.")


def main(
    phases: Annotated[Path,typer.Argument(exists=True,dir_okay=False)],
    native: Annotated[Path | None,typer.Option(exists=True,dir_okay=False)] = None,
):
    """Render measured phase and optional native tables to stdout."""
    data = json.loads(phases.read_text())
    for mode in ("nominal","full"):
        phase_table(mode,data[mode])
    if native:
        native_tables(json.loads(native.read_text()))


if __name__ == "__main__":
    typer.run(main)
