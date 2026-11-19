"""Release-free independent NNJVT preparation checks on original PHYSLITE."""
from __future__ import annotations

import ctypes
import hashlib
import json
import os
from pathlib import Path
import sys

import awkward as ak
import numpy as np
import typer

from histogram_reference import root_file

app = typer.Typer(add_completion=False)
JETS = "AnalysisJetsAuxDyn."
VERTICES = "PrimaryVerticesAuxDyn.vertexType"
FIELDS = [JETS + x for x in ("pt", "eta", "SumPtTrkPt500", "JVFCorr")] + [VERTICES]
FLT = np.ctypeslib.ndpointer(dtype=np.float64, flags=("C_CONTIGUOUS",))
INT = np.ctypeslib.ndpointer(dtype=np.int64, flags=("C_CONTIGUOUS",))


def load_library(library: Path):
    native = ctypes.CDLL(str(library.absolute()))
    native.met_jvt_error.restype = ctypes.c_char_p
    native.met_jvt_open.restype = ctypes.c_void_p
    native.met_jvt_open.argtypes = [ctypes.c_char_p]
    native.met_jvt_close.argtypes = [ctypes.c_void_p]
    native.met_jvt_evaluate.argtypes = [ctypes.c_void_p, ctypes.c_size_t, FLT, FLT, FLT, FLT, INT, INT, FLT]
    native.met_jvt_evaluate.restype = ctypes.c_int
    return native


def features(vertices, tracks) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    track_input, vertex_kind, feature_kind = [], [], []
    for types, jets in zip(ak.to_list(vertices), ak.to_list(tracks), strict=True):
        if 1 in types:
            pv = types.index(1)
            for moments in jets:
                if pv >= len(moments):
                    raise ValueError("NNJVT update has no track sum for the first primary vertex.")
                track_input.append(moments[pv])
                vertex_kind.append(2)
                feature_kind.append(1)
        elif types == [0]:
            track_input.extend([0.0] * len(jets))
            vertex_kind.extend([1] * len(jets))
            feature_kind.extend([0] * len(jets))
        else:
            raise ValueError(f"NNJVT has no primary or single dummy vertex in {types}.")
    return np.asarray(track_input, dtype=np.float64), np.asarray(vertex_kind, dtype=np.int64), np.asarray(feature_kind, dtype=np.int64)


def compute(native, handle, chunk):
    pt = ak.to_numpy(ak.flatten(chunk[JETS + "pt"])).astype(np.float64)
    eta = ak.to_numpy(ak.flatten(chunk[JETS + "eta"])).astype(np.float64)
    jvf = ak.to_numpy(ak.flatten(chunk[JETS + "JVFCorr"])).astype(np.float64)
    sums, vertex, available = features(chunk[VERTICES], chunk[JETS + "SumPtTrkPt500"])
    assert len(pt) == len(sums) == len(jvf) == len(eta)
    output = np.empty((len(pt), 3), dtype=np.float64)
    if native.met_jvt_evaluate(handle, len(pt), pt, eta, sums, jvf, vertex, available, output):
        raise RuntimeError(native.met_jvt_error().decode())
    return ak.to_numpy(ak.num(chunk[JETS + "pt"])), output


@app.command()
def calculate(
    input_file: Path, library: Path, payload: Path, result_path: Path,
    events: int = typer.Option(64, min=1, max=20000),
) -> None:
    if sys.version_info[:2] != (3, 14) or os.environ.get("AtlasVersion"):
        raise RuntimeError("Use a root-free Python 3.14 environment.")
    result_path.mkdir()
    native = load_library(library)
    handle = native.met_jvt_open(os.fsencode(payload))
    if not handle:
        raise RuntimeError(native.met_jvt_error().decode())
    counts, result = [], []
    try:
        with root_file(input_file) as original:
            source = original["CollectionTree"]
            for chunk in source.iterate(FIELDS, entry_stop=events, step_size=2000, library="ak"):
                n, values = compute(native, handle, chunk)
                counts.append(n)
                result.append(values)
            uuid = str(original.file.uuid)
    finally:
        native.met_jvt_close(handle)
    counts = np.concatenate(counts)
    values = np.concatenate(result)
    offsets = np.concatenate(([0], np.cumsum(counts)))
    assert offsets[-1] == len(values) and len(counts) == events
    with (result_path / "nnjvt.npz").open("xb") as output:
        np.savez(output, event_offsets=offsets, score=values[:, 0], decorator_pass=values[:, 1], met_accept=values[:, 2])
    mapped = sorted({line.split(maxsplit=5)[-1] for line in Path("/proc/self/maps").read_text().splitlines() if "correctionlib.so" in line or "libmet_nnjvt_probe.so" in line})
    maps = Path("/proc/self/maps").read_text()
    assert "libCore.so" not in maps and "/AnalysisBase/" not in maps
    report = {
        "source_fields": FIELDS, "input_uuid": uuid, "events": events, "jets": len(values),
        "python": sys.version, "process_threads": len(list(Path("/proc/self/task").iterdir())),
        "cpu_affinity": sorted(os.sched_getaffinity(0)), "loaded_probe_and_core": mapped,
        "payload_sha256": hashlib.sha256(payload.read_bytes()).hexdigest(),
    }
    with (result_path / "runtime.json").open("x") as output:
        json.dump(report, output, indent=2)


def stock_selector(trace: Path, events: int) -> np.ndarray:
    decisions = []
    with trace.open() as file:
        for row in file:
            entry, state, action, *fields = row.rstrip("\n").split("\t")
            if int(entry) >= events:
                break
            if state == "0" and action == "jetJvt":
                decisions.append(int(fields[-1]))
    return np.asarray(decisions)


@app.command()
def compare(result_path: Path, stock_ntuple: Path, stock_trace: Path, output_file: Path) -> None:
    with np.load(result_path / "nnjvt.npz") as actual, root_file(stock_ntuple) as stock:
        obs = stock["analysis"]
        events = len(actual["event_offsets"]) - 1
        score = obs["jet_NNJvt"].array(entry_stop=events)
        flags = obs["jet_NNJvtPass"].array(entry_stop=events)
        np.testing.assert_array_equal(ak.to_numpy(ak.num(score)), np.diff(actual["event_offsets"]))
        expected_score = ak.to_numpy(ak.flatten(score))
        expected_flags = ak.to_numpy(ak.flatten(flags))
        expected_accept = stock_selector(stock_trace, events)
        score_nonmatching = np.flatnonzero(expected_score != actual["score"])
        decorator_flag_nonmatching = int(np.count_nonzero(expected_flags != actual["decorator_pass"]))
        met_flag_nonmatching = int(np.count_nonzero(expected_accept != actual["met_accept"]))
        report = {
            "events": events, "jets": len(expected_score),
            "score_nonmatching": len(score_nonmatching),
            "score_max_abs": float(np.max(abs(expected_score - actual["score"]), initial=0)),
            "first_score_nonmatch": score_nonmatching[:10].tolist(),
            "decorator_flag_nonmatching": decorator_flag_nonmatching,
            "met_flag_nonmatching": met_flag_nonmatching,
        }
    with output_file.open("x") as output:
        json.dump(report, output, indent=2)
        output.write("\n")


if __name__ == "__main__":
    app()
