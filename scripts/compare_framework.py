"""Compare a real B capture directly with independent stock A."""
from __future__ import annotations

import json
from pathlib import Path

import awkward as ak
import numpy as np
import typer

from histogram_reference import met_gev, root_file
from summarize_observer import FINAL_QUANTITIES

app = typer.Typer(add_completion=False)
COMPONENTS = ("met", "mpx", "mpy", "sumet")


def numerical_quantity(branch: str, states: list[str]) -> tuple[str, str] | None:
    for quantity in (*COMPONENTS, "phi", *FINAL_QUANTITIES):
        for state in states:
            if branch in (f"finalMet_{quantity}_{state}", f"fullMet_{quantity}_{state}"):
                return quantity, state
    return None


def compare_numeric(original, candidate, quantity: str) -> dict:
    expected_storage = ak.to_numpy(original)
    actual_storage = ak.to_numpy(candidate)
    if expected_storage.shape != actual_storage.shape or expected_storage.dtype != actual_storage.dtype:
        raise ValueError(
            f"{quantity} storage differs: {expected_storage.shape} {expected_storage.dtype} "
            f"versus {actual_storage.shape} {actual_storage.dtype}"
        )
    expected = expected_storage.astype(np.float64)
    actual = actual_storage.astype(np.float64)
    exact = expected == actual
    equal_classification = np.isfinite(expected) == np.isfinite(actual)
    for predicate in (np.isnan, np.isposinf, np.isneginf):
        equal_classification &= predicate(expected) == predicate(actual)
    finite = np.isfinite(expected) & np.isfinite(actual)
    delta = np.zeros_like(expected)
    delta[finite] = actual[finite] - expected[finite]
    if quantity == "phi":
        delta[finite] = np.remainder(delta[finite] + np.pi, 2 * np.pi) - np.pi
        local_atol, local_rtol = 5e-7, 0.0
    elif quantity in COMPONENTS:
        local_atol, local_rtol = 0.001, 5e-7
    else:
        local_atol, local_rtol = 1e-6, 5e-7
    violation = ~equal_classification
    violation[finite] |= abs(delta[finite]) > (local_atol + local_rtol * abs(expected[finite]))
    failed = np.argwhere(violation)
    return {
        "shape": list(expected_storage.shape), "storage_dtype": str(expected_storage.dtype),
        "values": expected.size, "exact_finite": int(np.count_nonzero(exact & finite)),
        "max_absolute_difference": float(np.max(abs(delta), initial=0)),
        "nonfinite_classification_failures": int(np.count_nonzero(~equal_classification)),
        "total_failures": len(failed), "first_failure_positions": failed[:8].tolist(),
    }


@app.command()
def compare(
    stock_dir: Path, candidate_dir: Path, output_file: Path,
    expected_events: int | None = None,
) -> None:
    with root_file(stock_dir / "hist-reference.root") as stock_h, root_file(candidate_dir / "hist-candidate.root") as cand_h:
        states = stock_h["systematicsList"].axis().labels()
        assert states == cand_h["systematicsList"].axis().labels()
    with root_file(stock_dir / "reference.root") as stock_f, root_file(candidate_dir / "candidate.root") as cand_f:
        stock, candidate = stock_f["analysis"], cand_f["analysis"]
        n_events = candidate.num_entries
        stock_fields = stock.keys(recursive=False)
        cand_fields = candidate.keys(recursive=False)
        assert set(stock_fields) == set(cand_fields)
        assert 0 < n_events <= stock.num_entries
        if expected_events is not None:
            assert n_events == expected_events, (n_events, expected_events)
        identifiers = stock["eventNumber"].array(entry_stop=n_events, library="np") if "eventNumber" in stock else None
        report = {
            "events": n_events, "state_order": states,
            "numerical_checks": {}, "discrete_failures": [], "histogram_per_event_bin_failures": {},
        }
        for branch in stock_fields:
            if stock[branch].typename != candidate[branch].typename:
                raise ValueError(f"ROOT storage types differ for {branch}")
            expected = stock[branch].array(entry_stop=n_events)
            given = candidate[branch].array()
            key = numerical_quantity(branch, states)
            if key:
                report["numerical_checks"][branch] = compare_numeric(expected, given, key[0])
            elif not ak.almost_equal(expected, given, rtol=0, atol=0):
                report["discrete_failures"].append(branch)
        edges = np.linspace(0, 500, 101)
        for state in states:
            branch = f"finalMet_met_{state}"
            expected = met_gev(stock[branch].array(entry_stop=n_events, library="np"))
            observed = met_gev(candidate[branch].array(library="np"))
            changed = np.flatnonzero(np.searchsorted(edges, expected, side="right") != np.searchsorted(edges, observed, side="right"))
            report["histogram_per_event_bin_failures"][state] = {
                "count": len(changed), "first_failure_entries": changed[:8].tolist(),
                "first_failure_event_numbers": identifiers[changed[:8]].tolist() if identifiers is not None else [],
            }
    report["passed"] = (
        not report["discrete_failures"]
        and all(entry["total_failures"] == 0 for entry in report["numerical_checks"].values())
        and all(entry["count"] == 0 for entry in report["histogram_per_event_bin_failures"].values())
    )
    with output_file.open("x") as output:
        json.dump(report, output, indent=2, allow_nan=False)
        output.write("\n")
    if not report["passed"]:
        raise typer.Exit(1)


if __name__ == "__main__":
    app()
