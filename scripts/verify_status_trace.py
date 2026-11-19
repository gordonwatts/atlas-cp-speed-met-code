"""Check source-level stock statuses against the uninstrumented stock result."""
from __future__ import annotations

from collections import Counter
import ctypes
from itertools import groupby
import json
from pathlib import Path

import awkward as ak
import numpy as np
import typer
import uproot

from summarize_observer import TERMS

app = typer.Typer(add_completion=False)
ALGORITHM_ACTIONS = ("makerAlg", "builderAlg", "significanceAlg")
CODE_ACTIONS = ("rebuildMET", "rebuildJetMET", "applySystematicVariation", "applyCorrection", "varianceMET")
HARD_PARTICLES = ("RefEle", "RefGamma", "RefTau", "Muons")


def tree_open(result_path: Path):
    return uproot.open(result_path / "reference.root", num_workers=1, num_fallback_workers=1)


def check_stock_parity(original, traced) -> int:
    assert set(traced.keys()) == set(original.keys())
    count = 0
    for key in traced.keys(recursive=False):
        reference = original[key].array(entry_stop=traced.num_entries)
        actual = traced[key].array()
        assert reference.type == actual.type, key
        if not ak.almost_equal(reference, actual, rtol=0, atol=0):
            raise AssertionError(f"Instrumentation changed stock branch {key}")
        count += 1
    return count


def verify_rows(
    rows: list, event: int, states: list, soft_by_state: dict, soft_phi_nominal, jet_counts, summary: dict
) -> None:
    per_state = [dict() for _ in states]
    for state_index, action, fields in rows:
        if action in ALGORITHM_ACTIONS:
            assert state_index == len(states)
            assert fields == ["1"], (event, action, fields)
            summary["algorithm_calls"][action] += 1
            continue
        assert 0 <= state_index < len(states)
        per_state[state_index].setdefault(action, []).append(fields)
    baseline_seed = None
    baseline_draws = []
    baseline_inputs = None
    for state_index, (label, actions) in enumerate(zip(states, per_state, strict=True)):
        report = summary["states"][label]
        for action, calls in actions.items():
            report["calls_by_action"][action] += len(calls)
        assert [row[0] for row in actions["rebuildMET"]] == list(HARD_PARTICLES)
        assert all(row[-1] == "1" for row in actions["rebuildMET"])
        for term, count, _ in actions["rebuildMET"]:
            report["selected_leptons_by_term"][term] += int(count)
        assert actions["rebuildJetMET"] == [["RefJet", "PVSoftTrk", str(jet_counts[event]), "1", "1"]]
        assert actions["applySystematicVariation"] == [[label if label != "NOSYS" else "", "1"]]
        assert actions["applyCorrection"][0][0] == "2", (event, label, actions["applyCorrection"])
        assert actions["varianceMET"][0][0] == "1", (event, label, actions["varianceMET"])
        assert len(actions["varianceMET"]) == len(actions["applyCorrection"]) == 1
        soft = list(map(float.fromhex, actions["applyCorrection"][0][1:]))
        expected = [float(arr[event]) for arr in soft_by_state[label]]
        assert soft[3:] == expected, (event, label, soft[3:], expected)
        seed = int(actions["setRandomSeed"][0][0])
        expected_seed = ctypes.c_ulong(int(1e6 * float(soft_phi_nominal[event]))).value
        assert seed == expected_seed, (event, label, seed, expected_seed)
        root_seed_calls = actions["rootSeed"]
        expected_root_calls = [[str(seed or 1)]]
        if event == 0 and state_index == 0:
            # The lazily constructed TRandom3 seeds itself before the first set.
            assert root_seed_calls[0] == ["4357"]
            report["lazy_constructor_seed"] = 4357
            expected_root_calls.insert(0, ["4357"])
        assert root_seed_calls == expected_root_calls, (event, label, root_seed_calls)
        report["zero_seed_fallbacks"] += seed == 0
        report["negative_signed_phi_seeds"] += seed > 2 ** (ctypes.sizeof(ctypes.c_ulong) * 8 - 1)
        draws = actions.get("gaussian", [])
        if label == "NOSYS":
            baseline_seed = seed
            assert not draws
            assert soft[:3] == soft[3:]
        else:
            assert seed == baseline_seed
            assert len(draws) in (1, 2)
            assert draws[0][:2] == ["0x0p+0", "0x1p+0"]
            if baseline_draws:
                assert draws[0] == baseline_draws[0]
            else:
                baseline_draws = draws
        for _, local_index, passed in actions.get("associationFlag", []):
            assert passed in ("0", "1")
            local_index = int(local_index)
            current = summary["maximum_executed_association_flag_slot"]
            summary["maximum_executed_association_flag_slot"] = max(current, local_index)
        jvt = actions["jetJvt"] if jet_counts[event] else actions.get("jetJvt", [])
        assert len(jvt) == jet_counts[event], (event, label, len(jvt), jet_counts[event])
        assert [int(row[0]) for row in jvt] == list(range(jet_counts[event]))
        report["jet_jvt_accept"] += sum(row[1] == "1" for row in jvt)
        inputs = (soft[:3], jvt, actions.get("associationFlag", []), actions["rebuildMET"])
        if label == "NOSYS":
            baseline_inputs = inputs
        else:
            assert inputs == baseline_inputs, (event, label)


@app.command()
def verify(result_path: Path, stock_path: Path, output_path: Path) -> None:
    with uproot.open(result_path / "hist-reference.root") as hist:
        states = hist["systematicsList"].axis().labels()
    with tree_open(stock_path) as stock, tree_open(result_path) as traced:
        original, trace_tree = stock["analysis"], traced["analysis"]
        result = {
            "entries": trace_tree.num_entries,
            "states": {
                state: {
                    "calls_by_action": Counter(), "zero_seed_fallbacks": 0,
                    "negative_signed_phi_seeds": 0, "jet_jvt_accept": 0,
                    "selected_leptons_by_term": Counter(),
                }
                for state in states
            },
            "algorithm_calls": Counter(),
            "maximum_executed_association_flag_slot": -1,
            "physics_branch_payloads_exact_to_stock": check_stock_parity(original, trace_tree),
        }
        soft_index = list(TERMS).index("PVSoftTrk")
        soft_by_state = {}
        for state in states:
            soft_by_state[state] = [
                ak.to_numpy(trace_tree[f"fullMet_{component}_{state}"].array()[:, soft_index])
                for component in ("mpx", "mpy", "sumet")
            ]
        soft_phi_nominal = ak.to_numpy(trace_tree["fullMet_phi_NOSYS"].array()[:, soft_index])
        jet_counts = ak.to_list(ak.num(trace_tree["jet_pt"].array()))
        with (result_path / "status.tsv").open() as trace:
            lines = (line.rstrip("\n").split("\t") for line in trace)
            observed_entries = 0
            for entry_text, event_rows in groupby(lines, key=lambda row: row[0]):
                assert int(entry_text) == observed_entries
                rows = [(int(row[1]), row[2], row[3:]) for row in event_rows]
                verify_rows(rows, observed_entries, states, soft_by_state, soft_phi_nominal, jet_counts, result)
                observed_entries += 1
        assert observed_entries == trace_tree.num_entries
        assert dict(result["algorithm_calls"]) == dict.fromkeys(ALGORITHM_ACTIONS, observed_entries)
        result["trace_rows_aligned_with_entries_and_states"] = True
        result["status_code"] = {"SUCCESS": 1, "FAILURE": 0, "RECOVERABLE": 2}
        result["correction_code"] = {"Ok": 2, "Error": 0, "OutOfValidityRange": 1}
        result["remaining_unobservables"] = [
            "Individual private significance subtool return codes and covariance intermediates",
            "Uncalled options and deliberately unexercised malformed-input failure paths",
        ]
    with output_path.open("x") as output:
        json.dump(result, output, indent=2)
        output.write("\n")


if __name__ == "__main__":
    app()
