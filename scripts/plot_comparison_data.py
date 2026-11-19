"""Save compact matched distributions and residuals for every MET observable."""
from __future__ import annotations

import json
from pathlib import Path

import awkward as ak
import numpy as np
import typer

from histogram_reference import root_file
from summarize_observer import TERMS

app = typer.Typer(add_completion=False)
AXES = {
    "mpx": (-500000., 500000., -25000., 25000., "MeV"),
    "mpy": (-500000., 500000., -25000., 25000., "MeV"),
    "sumet": (0., 2000000., -25000., 25000., "MeV"),
    "met": (0., 500000., -25000., 25000., "MeV"),
    "phi": (-np.pi, np.pi, -1., 1., "rad"),
    "significance": (0., 40., -3., 3., "dimensionless"),
    "sigDirectional": (0., 40., -3., 3., "dimensionless"),
    "METOverSqrtSumET": (0., 25., -3., 3., "sqrt(GeV)"),
    "METOverSqrtHT": (0., 25., -3., 3., "sqrt(GeV)"),
}
N_BINS = 120


def wrapped_delta(current, nominal, is_phi: bool):
    delta = current - nominal
    if is_phi:
        delta = np.remainder(delta + np.pi, 2 * np.pi) - np.pi
    return delta


def histogram(array, edges) -> np.ndarray:
    if array.ndim == 1:
        array = array[:, np.newaxis]
    columns = []
    for values in array.T:
        values = values[np.isfinite(values)]
        indices = np.searchsorted(edges, values, side="right")
        columns.append(np.bincount(indices, minlength=len(edges) + 1))
    return np.stack(columns)


def summarize(array) -> dict:
    if array.ndim == 1:
        array = array[:, np.newaxis]
    rows = []
    for values in array.T:
        finite = abs(values[np.isfinite(values)])
        rows.append({
            "finite": len(finite), "nonfinite": len(values) - len(finite),
            "exact": int(np.count_nonzero(finite == 0)),
            "absolute_quantiles": np.quantile(finite, [0.5, 0.9, 0.99]).tolist() if len(finite) else [None] * 3,
            "max_absolute": float(np.max(finite, initial=0)),
        })
    return rows


@app.command()
def prepare(
    a_file: Path, candidate_file: Path, candidate_metadata: Path,
    output_npz: Path, output_json: Path,
) -> None:
    with root_file(candidate_metadata) as meta:
        states = meta["systematicsList"].axis().labels()
    with root_file(a_file) as stock, root_file(candidate_file) as trial:
        a, b = stock["analysis"], trial["analysis"]
        n = b.num_entries
        report = {
            "events": n, "terms": list(TERMS), "states": states,
            "fields": {}, "residual_quantiles": [0.5, 0.9, 0.99], "density_flow_bins": True,
        }
        plot = {"states": np.asarray(states), "terms": np.asarray(list(TERMS))}
        for quantity, (low, high, delta_low, delta_high, units) in AXES.items():
            prefix = "fullMet" if quantity in list(AXES)[:5] else "finalMet"
            edges = np.linspace(low, high, N_BINS + 1)
            delta_edges = np.linspace(delta_low, delta_high, N_BINS + 1)
            plot[f"{quantity}_edges"] = edges
            plot[f"{quantity}_delta_edges"] = delta_edges
            keys = [f"{prefix}_{quantity}_{state}" for state in states]
            expected = [ak.to_numpy(a[key].array(entry_stop=n)).astype(np.float64) for key in keys]
            actual = [ak.to_numpy(b[key].array()).astype(np.float64) for key in keys]
            phi = quantity == "phi"
            report["fields"][quantity] = {
                "units": units, "terms": list(TERMS) if prefix == "fullMet" else ["Final"],
                "comparison": {},
            }
            for source, rows in (("a", expected), ("candidate", actual)):
                plot[f"{quantity}_{source}"] = np.stack([histogram(values, edges) for values in rows])
                plot[f"{quantity}_{source}_delta"] = np.stack([
                    histogram(wrapped_delta(values, rows[0], phi), delta_edges) for values in rows
                ])
            for label, a_state, b_state in zip(states, expected, actual, strict=True):
                residual = wrapped_delta(b_state, a_state, phi)
                a_shift = wrapped_delta(a_state, expected[0], phi)
                b_shift = wrapped_delta(b_state, actual[0], phi)
                report["fields"][quantity]["comparison"][label] = {
                    "residual": summarize(residual),
                    "delta_from_nominal_residual": summarize(wrapped_delta(b_shift, a_shift, phi)),
                }
    with output_npz.open("xb") as output:
        np.savez_compressed(output, **plot)
    with output_json.open("x") as output:
        json.dump(report, output, indent=2, allow_nan=False)
        output.write("\n")


if __name__ == "__main__":
    app()
