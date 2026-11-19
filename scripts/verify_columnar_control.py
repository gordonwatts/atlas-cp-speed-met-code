"""Compare independently recorded C entry/bin assignments with the stock A ntuple."""
from pathlib import Path
import json

import numpy as np
import typer

from histogram_reference import met_gev, root_file

app = typer.Typer(add_completion=False)


@app.command()
def verify(
    stock_result: Path, columnar_result: Path, output_file: Path,
) -> None:
    a_hist = stock_result / "histogram.npz"
    c_hist = columnar_result / "histogram.npz"
    assert a_hist.read_bytes() == c_hist.read_bytes()
    with np.load(a_hist) as a_saved, np.load(columnar_result / "entry-bins.npz") as c_saved:
        edges = a_saved["bin_edges_gev"]
        with root_file(stock_result / "reference.root") as stock:
            arrays = stock["analysis"].arrays(["core_mpx", "core_mpy"], library="np")
        gev = np.hypot(met_gev(arrays["core_mpx"]), met_gev(arrays["core_mpy"]))
        a_bins = np.searchsorted(edges, gev, side="right")
        np.testing.assert_array_equal(np.arange(len(a_bins)), c_saved["entry"])
        np.testing.assert_array_equal(a_bins, c_saved["flow_inclusive_bin"])
    report = {
        "observable": "PVSoftTrkCore",
        "original_fields_C": ["name", "mpx", "mpy"],
        "entries": len(a_bins),
        "flow_inclusive_ordered_entry_bin_assignments_match": True,
        "serialized_npz_equal": True,
    }
    with output_file.open("x") as output:
        json.dump(report, output, indent=2)
        output.write("\n")


if __name__ == "__main__":
    app()
