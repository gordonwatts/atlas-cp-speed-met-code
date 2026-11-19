"""Independently compare standalone MET float storage and output bins with A."""
from __future__ import annotations

import json
from pathlib import Path

import awkward as ak
import numpy as np
import typer

from compare_framework import compare_numeric
from histogram_reference import met_gev, root_file
from plot_comparison_data import AXES, N_BINS, histogram, summarize, wrapped_delta
from summarize_observer import FINAL_QUANTITIES, TERMS

app = typer.Typer(add_completion=False)
TERM_FIELDS = ("mpx", "mpy", "sumet", "met", "phi")


def plot_fields(expected: dict, actual: dict, states: list, plot_path: Path, report_path: Path) -> None:
    plot = {"states": np.asarray(states), "terms": np.asarray(list(TERMS))}
    report = {"fields": {}, "terms": list(TERMS), "states": states, "density_flow_bins": True}
    for field, (low, high, delta_low, delta_high, units) in AXES.items():
        edges = np.linspace(low, high, N_BINS + 1)
        delta_edges = np.linspace(delta_low, delta_high, N_BINS + 1)
        plot[f"{field}_edges"] = edges
        plot[f"{field}_delta_edges"] = delta_edges
        is_phi = field == "phi"
        a, c = (np.asarray(group[field], dtype=np.float64) for group in (expected, actual))
        field_report = {"units": units, "terms": list(TERMS) if field in TERM_FIELDS else ["Final"], "comparison": {}}
        for name, rows in (("a", a), ("candidate", c)):
            plot[f"{field}_{name}"] = np.stack([histogram(row, edges) for row in rows])
            plot[f"{field}_{name}_delta"] = np.stack([histogram(wrapped_delta(row, rows[0], is_phi), delta_edges) for row in rows])
        for state, a_row, c_row in zip(states, a, c, strict=True):
            field_report["comparison"][state] = {
                "residual": summarize(wrapped_delta(c_row, a_row, is_phi)),
                "delta_from_nominal_residual": summarize(wrapped_delta(
                    wrapped_delta(c_row, c[0], is_phi), wrapped_delta(a_row, a[0], is_phi), is_phi
                )),
            }
        report["fields"][field] = field_report
    with plot_path.open("xb") as output:
        np.savez_compressed(output, **plot)
    with report_path.open("x") as output:
        json.dump(report, output, indent=2, allow_nan=False)
        output.write("\n")


@app.command()
def compare(stock_dir: Path, candidate: Path, expected_events: int) -> None:
    with root_file(stock_dir / "hist-reference.root") as metadata:
        states = metadata["systematicsList"].axis().labels()
    expected = {}
    actual = {}
    with (
        root_file(stock_dir / "reference.root") as stock_file,
        np.load(candidate / "c-observables.npz", allow_pickle=False) as observed,
        np.load(candidate / "c.npz", allow_pickle=False) as saved_hist,
    ):
        stock = stock_file["analysis"]
        n = len(observed["event_number"])
        assert n == expected_events <= stock.num_entries
        assert observed["states"].tolist() == states
        assert observed["term_names"].tolist() == list(TERMS)
        assert observed["term_fields"].tolist() == list(TERM_FIELDS)
        quantities = observed["terms"]
        sources = observed["sources"]
        significance = observed["significance"]
        assert quantities.shape == (n, len(states), len(TERMS), len(TERM_FIELDS)) and quantities.dtype == np.dtype("f4")
        assert sources.shape == (n, len(states), len(TERMS)) and sources.dtype == np.dtype("u8")
        assert significance.shape == (n, len(states), len(FINAL_QUANTITIES)) and significance.dtype == np.dtype("f4")
        expected_event = stock["eventNumber"].array(entry_stop=n, library="np")
        assert observed["event_number"].dtype == expected_event.dtype
        np.testing.assert_array_equal(observed["event_number"], expected_event)
        assert observed["entry"].dtype == np.dtype("u8")
        np.testing.assert_array_equal(observed["entry"], np.arange(n, dtype="u8"))
        coverage = {"event_number": True, "entry": True}
        for field, stock_name in (("run_number", "runNumber"), ("mc_channel_number", "mcChannelNumber")):
            expected_identity = stock[stock_name].array(entry_stop=n, library="np")
            assert observed[field].dtype == expected_identity.dtype
            np.testing.assert_array_equal(observed[field], expected_identity)
            coverage[field] = True
        report = {
            "events": n, "states": states, "identity_observation_coverage": coverage,
            "states_checked": {}, "first_bin_event_mismatch": None,
            "numerical_validation_passed": True,
        }
        for quantity in (*TERM_FIELDS, *FINAL_QUANTITIES):
            expected[quantity] = []
            actual[quantity] = []
        edges = np.linspace(0, 500, 101)
        assert saved_hist["state_labels"].tolist() == states
        np.testing.assert_array_equal(saved_hist["bin_edges_gev"], edges)
        counts = []
        for index, state in enumerate(states):
            state_report = {}
            report["states_checked"][state] = state_report
            a_masks = ak.to_numpy(stock[f"fullMet_source_{state}"].array(entry_stop=n))
            assert a_masks.shape == sources[:, index, :].shape and a_masks.dtype == sources.dtype
            np.testing.assert_array_equal(a_masks, sources[:, index, :])
            a_names = stock[f"fullMet_term_{state}"].array(entry_stop=n)
            assert ak.to_list(a_names) == [list(TERMS)] * n
            for field_index, field in enumerate(TERM_FIELDS):
                expected[field].append(ak.to_numpy(stock[f"fullMet_{field}_{state}"].array(entry_stop=n)))
                actual[field].append(quantities[:, index, :, field_index])
                state_report[field] = compare_numeric(expected[field][-1], actual[field][-1], field)
            for field_index, field in enumerate(FINAL_QUANTITIES):
                expected[field].append(stock[f"finalMet_{field}_{state}"].array(entry_stop=n, library="np"))
                actual[field].append(significance[:, index, field_index])
                state_report[field] = compare_numeric(expected[field][-1], actual[field][-1], field)
            reference_magnitude = expected["met"][-1][:, -1]
            output_magnitude = actual["met"][-1][:, -1]
            assert np.isfinite(reference_magnitude).all()
            assert np.isfinite(output_magnitude).all()
            stock_bins = np.searchsorted(edges, met_gev(reference_magnitude), side="right")
            output_bins = np.searchsorted(edges, met_gev(output_magnitude), side="right")
            changed = np.flatnonzero(stock_bins != output_bins)
            if len(changed) and report["first_bin_event_mismatch"] is None:
                report["first_bin_event_mismatch"] = {
                    "state": state, "entry": int(changed[0]), "event_number": int(expected_event[changed[0]]),
                }
            if len(changed) or any(row["total_failures"] for row in state_report.values()):
                report["numerical_validation_passed"] = False
            state_report["histogram"] = {
                "per_event_bin_mismatches": len(changed), "stock_flow": [int(np.sum(stock_bins == 0)), int(np.sum(stock_bins == 101))],
                "candidate_flow": [int(np.sum(output_bins == 0)), int(np.sum(output_bins == 101))],
            }
            counts.append(np.bincount(stock_bins, minlength=102))
        histogram_equal = all(np.array_equal(saved_hist[key], np.stack(counts)) for key in ("counts", "sumw", "sumw2"))
        report["saved_histogram_equals_stock"] = histogram_equal
        report["numerical_validation_passed"] &= histogram_equal
    with (candidate / "comparison.json").open("x") as output:
        json.dump(report, output, indent=2, allow_nan=False)
        output.write("\n")
    plot_fields(expected, actual, states, candidate / "plot-comparison.npz", candidate / "plot-summary.json")
    if not report["numerical_validation_passed"]:
        raise typer.Exit(1)


if __name__ == "__main__":
    app()
