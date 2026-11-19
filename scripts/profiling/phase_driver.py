"""Instrument coarse scopes in an immutable copy of the accepted visible MET runner."""

import ast
from collections import defaultdict
from contextlib import contextmanager
import json
import os
from pathlib import Path
import sys
import time

METER = "_met_profile_recorder"
OBJECTS = ("jet", "electron", "photon", "tau", "muon")


def clock():
    return time.monotonic_ns(), time.process_time_ns()


def phase(label, body):
    result = ast.With(items=[ast.withitem(ast.Call(
        ast.Attribute(ast.Name(METER, ast.Load()), "phase", ast.Load()),
        [ast.Constant(label)], [],
    ))], body=body)
    return ast.copy_location(result, body[0])


def node(stmts, predicate, description):
    matches = [stmt for stmt in stmts if predicate(stmt)]
    if len(matches) != 1:
        raise ValueError(f"Expected exactly one {description}, found {len(matches)}")
    return matches[0]


def prefix(stmts, text):
    return node(stmts, lambda stmt: ast.unparse(stmt).startswith(text), text)


def replace(stmts, stmt, label):
    stmts[stmts.index(stmt)] = phase(label, [stmt])


def payload_call(call, label):
    return ast.copy_location(ast.Call(
        ast.Attribute(ast.Name(METER, ast.Load()), "call", ast.Load()),
        [label, call.func, *call.args], call.keywords,
    ), call)


def instrument(source, detailed):
    tree = ast.parse(source.read_bytes(), filename=str(source))
    main = node(tree.body, lambda item: isinstance(item, ast.FunctionDef) and item.name == "main", "visible main")
    outer = node(main.body, lambda item: isinstance(item, ast.Try) and bool(item.finalbody), "native context try/finally")
    processing = prefix(outer.body, "for _ in range(repeat):")
    loops = processing.body
    input_scope = prefix(loops, "with uproot.open(input_file")
    booking = prefix(loops, "histogram = book(labels)")
    digest = prefix(loops, "digest = signature(histogram)")
    digest_result = prefix(loops, "expected_digest = digest")
    if loops[-1] is not digest_result or not isinstance(input_scope, ast.With):
        raise ValueError("MET runner per-pass histogram boundary changed")
    chunks = prefix(input_scope.body, "for first in range(0, end, chunk_size):")
    packing = prefix(chunks.body, "buffers = pack_raw(raw, assoc, link_reader.read(")
    packing_index = chunks.body.index(packing)
    mark = ast.Expr(ast.Call(
        ast.Attribute(ast.Name(METER, ast.Load()), "pass_result", ast.Load()),
        [ast.Name("entries", ast.Load()), ast.Name("digest", ast.Load())], [],
    ))
    if detailed:
        # Keep the accepted script's original __future__ imports and insert a
        # top-level timing scope around its standard and third-party imports.
        begin = max(i for i, stmt in enumerate(tree.body)
                    if isinstance(stmt, ast.ImportFrom) and stmt.module == "__future__") + 1
        uuid = tree.body.index(prefix(tree.body, "UUID = "))
        tree.body[begin:uuid] = [phase("imports", tree.body[begin:uuid])]
        outer_index = main.body.index(outer)
        main.body[:outer_index] = [phase("initialize", main.body[:outer_index])]
        replace(outer.body, prefix(outer.body, "bridge.configure("), "initialize")
        replace(outer.finalbody, prefix(outer.finalbody, "bridge.lib.met_columnar_close("), "close_native")
        replace(main.body, prefix(main.body, "save(histogram, labels, entries, output)"), "save_histogram")
        replace(loops, prefix(loops, "bridge.lib.met_columnar_reset("), "pass_reset")
        replace(loops, booking, "book_histogram")
        replace(loops, digest, "pass_histogram_signature")

        opener = input_scope.items[0].context_expr
        if len(input_scope.items) != 1 or not isinstance(opener, ast.Call):
            raise ValueError("Changed original file opener structure")
        input_scope.items[0].context_expr = payload_call(opener, ast.Constant("open_original"))
        header = input_scope.body.index(chunks)
        input_scope.body[:header] = [phase("input_metadata", input_scope.body[:header])]
        replace(chunks.body, prefix(chunks.body, "raw = tree.arrays("), "input_read_decode")
        replace(chunks.body, prefix(chunks.body, "assoc = read_associations("), "input_read_decode")
        jet_read = packing.value.args[2]
        packing.value.args[2] = payload_call(jet_read, ast.Constant("input_read_decode"))
        replace(chunks.body, packing, "raw_transport")

        evaluator = prefix(chunks.body, "evaluate = lambda name, inputs: bridge.prepare(")
        direct = ast.BinOp(ast.Constant("nominal_node."), ast.Add(), ast.Name("name", ast.Load()))
        evaluator.value.body = payload_call(evaluator.value.body, direct)
        replace(chunks.body, prefix(chunks.body, "prepare_chunk("), "nominal_projection")
        shared = node(chunks.body, lambda stmt: isinstance(stmt, ast.Assign) and
                      isinstance(stmt.value, ast.Call) and ast.unparse(stmt.value.func) == "bridge.full",
                      "shared native event evaluation")
        histogram = node(chunks.body, lambda stmt: isinstance(stmt, ast.For) and
                         ast.unparse(stmt.iter) == "enumerate(labels)", "histogram state loop")
        replace(chunks.body, shared, "shared_event_evaluation")
        replace(chunks.body, histogram, "fill_histogram")
    objects = ast.Expr(ast.Call(
        ast.Attribute(ast.Name(METER, ast.Load()), "chunk_objects", ast.Load()),
        [ast.Name("buffers", ast.Load())], [],
    ))
    chunks.body.insert(packing_index+1, ast.copy_location(objects, packing))
    loops.append(ast.copy_location(mark, digest_result))
    processing.body = [phase("processing_pass", loops)]
    return compile(ast.fix_missing_locations(tree), str(source), "exec")


class Measure:
    def __init__(self):
        self.stack = []
        self.passes = []
        self.objects = {owner:0 for owner in OBJECTS}
        self.exclusive = defaultdict(lambda: [0, 0, 0])
        self.inclusive = defaultdict(lambda: [0, 0, 0])
        self.worker_affinity = sorted(os.sched_getaffinity(0))

    @contextmanager
    def phase(self, label):
        frame = [0, 0]
        self.stack.append(frame)
        original_pass_count = len(self.passes)
        if label == "processing_pass":
            self.objects = dict.fromkeys(OBJECTS, 0)
        start = clock()
        try:
            yield
        finally:
            end = clock()
            duration = [b-a for a,b in zip(start, end)]
            self.stack.pop()
            if label == "processing_pass":
                if len(self.passes) > original_pass_count:
                    self.passes[-1].update(start_monotonic_ns=start[0], end_monotonic_ns=end[0],
                                           wall_ns=duration[0], cpu_ns=duration[1])
                # Pass windows are inclusive fences and excluded from the stage budget.
            else:
                for kind, span in ((self.inclusive[label], duration),
                                   (self.exclusive[label], [d-child for d,child in zip(duration, frame)])):
                    kind[0] += span[0]
                    kind[1] += span[1]
                    kind[2] += 1
            if self.stack:
                self.stack[-1][0] += duration[0]
                self.stack[-1][1] += duration[1]

    def call(self, label, function, /, *args, **kwargs):
        with self.phase(label):
            return function(*args, **kwargs)

    def chunk_objects(self, buffers):
        for owner in OBJECTS:
            self.objects[owner] += int(buffers[owner+".o"][-1])

    def pass_result(self, entries, digest):
        self.passes.append({"pass":len(self.passes)+1, "entries":entries,
                            "histogram_sha256":digest, "objects":self.objects.copy()})

    @staticmethod
    def rows(data):
        return {key:dict(zip(("wall_ns", "cpu_ns", "calls"), values))
                for key, values in sorted(data.items())}


def run(source, output, mode, arguments):
    """Internal subprocess callable. The original Typer parser receives its usual argv."""
    source, output = Path(source), Path(output)
    if mode not in ("passes", "phases"):
        raise ValueError("Unknown diagnostic measurement mode")
    compiled = instrument(source, mode == "phases")
    meter = Measure()
    sys.path.insert(0, str(source.parent))
    sys.argv = [str(source), *arguments]
    start = clock()
    success = False
    try:
        exec(compiled, {"__name__":"__main__", "__file__":str(source), METER:meter})
        success = True
    except SystemExit as error:
        success = error.code in (None, 0)
        raise
    finally:
        end = clock()
        span = {"wall_ns":end[0]-start[0], "cpu_ns":end[1]-start[1]}
        exclusive = meter.rows(meter.exclusive)
        residual = {k:span[k]-sum(row[k] for row in exclusive.values()) for k in span}
        payload = {"schema":1, "mode":mode, "success":success, "pid":os.getpid(),
                   "worker_affinity":meter.worker_affinity,
                   "argv":sys.argv, "script_span":span, "script_start_monotonic_ns":start[0],
                   "script_end_monotonic_ns":end[0], "worker_cpu_at_exec_finish_ns":end[1],
                   "exclusive":exclusive, "inclusive":meter.rows(meter.inclusive),
                   "unassigned_in_script":residual, "passes":meter.passes}
        with output.open("x", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, allow_nan=False)
            handle.write("\n")
