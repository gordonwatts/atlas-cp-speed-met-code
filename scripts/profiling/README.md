# Standalone MET phase and native diagnostics

Use the pip-installed interpreter and visible script from the final accepted portable C bundle. The child environment removes release paths and limits numerical libraries to one thread. Analysis runs on CPU 0. The controller and perf collector run on CPU 1 by default and record their separate affinity. The CPU0 topology receipt shows whether that choice shares a physical core. Run the diagnostics only after the matched production timings are complete.

The original-file UUID, 20,000 rows, state labels and ordered histogram arrays are checked for every child. The harness records the visible script SHA256, wheel SHA256, bundle receipt, original file identity, exact argv, observed worker affinity and threads, and mapped library hashes matched to the wheel. The runner verifies an identical histogram digest after every pass. The harness separately compares all histogram arrays to the accepted stock output.

Use absolute paths and fresh output directories. Set `stage` to the validated portable checkpoint and `nominal` to the accepted independently written one-state histogram. If the one-state histogram is missing, use the accepted stock nominal AnalysisBase run and the existing `scripts/histogram_reference.py` program to save it before starting the diagnostics.

```bash
cd /home/gwatts/atlas-cp-speed/met-20261119
stage="$PWD/results/columnar-m26-portable-001"  # use the actual accepted stage
chunk=5000  # use the measured stage's chunk size
py="$stage/pip-env/bin/python"
read -r input < config/input.txt
full="$PWD/results/reference-histogram-transport-001/histogram.npz"
nominal="$ACCEPTED_STOCK_NOMINAL_NPZ"
phase_result="$PWD/results/c-phase-diagnostic-001"
native_result="$PWD/results/c-native-diagnostic-001"
common=("$py" scripts/profiling/measure.py
        --python "$py" --demo "$stage/bundle/met_demo.py"
        --input "$input" --bundle-manifest "$stage/bundle/bundle.json"
        --analysis-cpu 0 --collector-cpu 1 --chunk-size "$chunk")
"${common[@]}" phases --out "$phase_result" --accepted-nominal "$nominal" --accepted-full "$full"
passes=$("${common[@]}" estimate "$phase_result/summary.json")
"${common[@]}" native --out "$native_result" --accepted "$full" --repeats "$passes" --pmu
"$py" scripts/profiling/render.py "$phase_result/summary.json" --native "$native_result/summary.json"
```

## Elapsed phase boundary

The phase command runs one direct subprocess and one diagnostic subprocess for each mode. The direct control launches the unchanged visible entry point with the specified original-entry chunk size. The diagnostic subprocess independently compiles an in-memory copy of the same script. Structural anchors fail if its accepted entry-point shape changed. The AST adds coarse phase contexts and per-pass receipts to the imports and original operations. Capture/trace/runtime-record options stay off. There is no patch to a production or installed source file.

- Linux `wait4` records elapsed and process tree user/system CPU. Direct elapsed includes interpreter startup, imports, initialization, original input and output close. In sampled mode the process tree includes the collector. The diagnostic Python process separately records its own CPU.
- Named phases cover imports, initialization and configure, original opening/metadata, reading and decoding associations and branch fields, raw packing, direct nominal projection and its named native calls, terms/significance evaluation, histogram book/fill/save, and native close. Nested scopes subtract their children from exclusive time. The phase receipt also publishes inclusive totals with the original overlap.
- The residual within the visible script includes context-manager boundaries, dispatch and loop bookkeeping, untimed optional checks that remain off, and timer work. The separate outside-script residual includes launcher startup/teardown, AST construction, Python diagnostic driver work, and phase-file writing. Sampled outside-script process-tree CPU includes perf. The direct control gives the whole real executable boundary.
- Each diagnostic processing pass includes native reset, histogram book, reopening PHYSLITE from entry zero, all chunks and histogram digest. Only the last histogram is saved. Read-only visits are counted from existing original typed-buffer offsets. Do not add input-output workload changes to increase a measurement.

## Native capture

The software event probe requires a numeric `cpu-clock:u`. The probes record returned value and evidence even when an unsupported PMU reports exit status zero. If the environment supports both hardware counters, `--pmu` adds a matched repeated real workload under `perf stat`. Retain raw counts and the workload receipt. Missing or unavailable events remain a capability observation.

Native capture runs three matched commands in fresh directories with the same `-n`. The first is direct, the second adds only per-pass boundaries and object visit counts, and the third records DWARF stacks with the same pass-only instrumentation. The overhead receipt compares direct to pass-only and pass-only to sampled. These are diagnostic comparisons. The production timing ledger carries its own uninstrumented repetitions.

The profile expects at least 60 seconds both of measured processing elapsed and of actual process CPU, excluding setup and final output. This conservative CPU floor also limits the risk of treating an I/O wait as demonstrated sampling time. If a workload falls short, retain the failed directory and use the suggested repeat count in a fresh directory. Reopening the original file warms cached pages and amortizes initialization. The retained raw profile covers the whole process. The published denominator selects only the recorded monotonic windows of the single processing worker on the expected CPU. Report discarded startup/teardown and other processes separately. Report period weights, samples, detected loss, resolved depths and unparsed lines. Inspect those records before treating symbol coverage as adequate.

## Evidence files

| Artifact | Contents |
| --- | --- |
| Phase `summary.json`, per-run `validated.json` and `phase_receipt.json` | Nanosecond exclusive/inclusive stage budgets, unassigned wall and CPU, stopwatch boundaries, per-pass counts and digests, saved outputs and matched direct overhead |
| Per-run `command.json`, `process.json` and `processes.jsonl` | Exact argv, wait4 receipt, actual worker and collector identities, sampled affinity/threads, mapped libraries |
| Native `perf_tool.json`, `counter_support.json`, `*.stat` and logs | Version, privilege and tested software/hardware support |
| Sampled `perf.data`, `physical_stacks.txt`, `build_ids.txt` | Original raw recording, decoded physical frames and mapped build IDs |
| Sampled `self.tsv`, `inclusive.tsv`, `stacks.tsv` | Whole-window denominator, disjoint leaf costs, overlapping inclusive symbol costs and call paths |
| Sampled `named_remainder_by_library.tsv`, `stage_by_stack.tsv`, native `summary.json` | Known omitted functions by executing library, explicit stack-recognized stage or unattributed, independent unresolved-leaf count, coverage and overhead |

Self cost assigns each CPU-clock observation to one executing leaf. Inclusive counts a symbol anywhere in the stack once per sample and overlaps across rows. The known omitted functions, top self functions and unresolved observations use the same denominator, and their periods sum exactly before percentage rounding. Native stage labels require stack anchors. Library identity alone cannot assign a Python/allocator sample to I/O or evaluation. Sample shares estimate sampled CPU cost and cannot be added to elapsed phase shares. `render.py` prints compact tables from recorded receipts for the supporting study.
