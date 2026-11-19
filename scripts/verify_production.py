"""Check the minimal A endpoint against the observer and count bins independently."""
from __future__ import annotations

import json
from pathlib import Path

import awkward as ak
import numpy as np
import typer

from histogram_reference import met_gev, root_file

app = typer.Typer(add_completion=False)


def values_by_state(input_file: Path, result_path: Path, observer_path: Path | None) -> dict:
    with root_file(result_path / "hist-reference.root") as meta:
        labels = meta["systematicsList"].axis().labels()
    with root_file(result_path / "reference.root") as production:
        tree = production["analysis"]
        if observer_path is None:
            assert labels == ["NOSYS"]
            assert set(tree.keys()) == {"core_mpx", "core_mpy"}
            with root_file(input_file) as raw_file:
                original = raw_file["CollectionTree"]
                names = original["MET_Core_AnalysisMETAuxDyn.name"].array(entry_stop=tree.num_entries)
                selected = names == "PVSoftTrkCore"
                assert ak.all(ak.sum(selected, axis=1) == 1)
                for component in ("mpx", "mpy"):
                    raw = original[f"MET_Core_AnalysisMETAuxDyn.{component}"].array(entry_stop=tree.num_entries)
                    observed = tree[f"core_{component}"].array(library="np")
                    np.testing.assert_array_equal(observed, ak.to_numpy(ak.firsts(raw[selected], axis=1)))
            return {
                "PVSoftTrkCore": np.hypot(
                    met_gev(tree["core_mpx"].array(library="np")), met_gev(tree["core_mpy"].array(library="np"))
                )
            }

        fields = [f"finalMet_met_{state}" for state in labels]
        assert set(tree.keys()) == set(fields)
        with root_file(observer_path / "hist-reference.root") as metadata:
            assert labels == metadata["systematicsList"].axis().labels()
        with root_file(observer_path / "reference.root") as observer:
            assert observer["analysis"].num_entries == tree.num_entries
            result = {}
            for label, field in zip(labels, fields, strict=True):
                actual = tree[field].array(library="np")
                expected = observer["analysis"][field].array(library="np")
                np.testing.assert_array_equal(actual, expected)
                result[label] = met_gev(actual)
            return result


def histogram_check(path: Path, values: dict) -> dict:
    with np.load(path) as data:
        expected_edges = np.linspace(0, 500, 101)
        assert data["state_labels"].tolist() == list(values)
        np.testing.assert_array_equal(data["bin_edges_gev"], expected_edges)
        rows = []
        for state, gev in values.items():
            assert np.isfinite(gev).all(), state
            assignments = np.searchsorted(expected_edges, gev, side="right")
            rows.append(np.bincount(assignments, minlength=len(expected_edges) + 1))
        for field in ("counts", "sumw", "sumw2"):
            np.testing.assert_array_equal(data[field], np.stack(rows))
        return {
            "state_labels": list(values),
            "bin_edges_GeV": expected_edges.tolist(),
            "bin_count": len(expected_edges) - 1,
            "weights": "unit",
            "per_state_counts": {key: row.tolist() for key, row in zip(values, rows, strict=True)},
            "per_state_totals": {key: int(row.sum()) for key, row in zip(values, rows, strict=True)},
            "per_state_flows": {key: [int(row[0]), int(row[-1])] for key, row in zip(values, rows, strict=True)},
            "histogram_counts_sumw_sumw2_exact": True,
        }


@app.command()
def verify(input_file: Path, result_path: Path, output_path: Path, observer_path: Path | None = None) -> None:
    values = values_by_state(input_file, result_path, observer_path)
    report = histogram_check(result_path / "histogram.npz", values)
    report["reference_equals_production_per_event"] = True
    with output_path.open("x") as output:
        json.dump(report, output, indent=2)
        output.write("\n")


if __name__ == "__main__":
    app()
