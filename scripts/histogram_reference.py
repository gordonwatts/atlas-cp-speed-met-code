"""Write the fixed unit-weight MET histogram from stock AnalysisBase output."""
from __future__ import annotations

from pathlib import Path

import hist
import numpy as np
import typer
import uproot

app = typer.Typer(add_completion=False)
BINS = 100
MIN_GEV = 0.0
MAX_GEV = 500.0


def root_file(path: Path):
    return uproot.open(path, handler=uproot.source.file.MemmapSource, num_workers=1, num_fallback_workers=1)


def met_gev(mev: np.ndarray) -> np.ndarray:
    return mev.astype(np.float64) / 1000.0


def book(labels: list[str]) -> hist.Hist:
    return hist.Hist(
        hist.axis.StrCategory(labels, overflow=False, name="state"),
        hist.axis.Regular(BINS, MIN_GEV, MAX_GEV, name="met_GeV", underflow=True, overflow=True),
        storage=hist.storage.Weight(),
    )


def fill(histogram: hist.Hist, state: str, gev: np.ndarray) -> None:
    if not np.isfinite(gev).all():
        raise ValueError(f"Nonfinite MET in {state}; record event identities and freeze policy first.")
    histogram.fill(state, gev, threads=1)


def save(histogram: hist.Hist, labels: list[str], entries: int, output_file: Path) -> None:
    sumw = histogram.values(flow=True)
    sumw2 = histogram.variances(flow=True)
    counts = sumw.astype(np.int64)
    np.testing.assert_array_equal(counts, sumw)
    np.testing.assert_array_equal(counts, sumw2)
    np.testing.assert_array_equal(counts.sum(axis=1), np.full(len(labels), entries))
    with output_file.open("xb") as output:
        np.savez(
            output, state_labels=np.array(labels), bin_edges_gev=histogram.axes[1].edges,
            counts=counts, sumw=sumw, sumw2=sumw2,
        )

@app.command()
def write(analysis_file: Path, metadata_file: Path, output_file: Path, stored: bool = False) -> None:
    with root_file(metadata_file) as meta:
        labels = meta["systematicsList"].axis().labels()
        if not labels or labels[0] != "NOSYS" or len(set(labels)) != len(labels):
            raise ValueError(f"Invalid systematic metadata: {labels}")
    if stored:
        if labels != ["NOSYS"]:
            raise ValueError(f"Unexpected stored MET control states: {labels}")
        labels = ["PVSoftTrkCore"]
    histogram = book(labels)
    fields = ["core_mpx", "core_mpy"] if stored else [f"finalMet_met_{state}" for state in labels]
    with root_file(analysis_file) as original:
        for chunk in original["analysis"].iterate(fields, step_size=2000, library="np"):
            if stored:
                fill(histogram, labels[0], np.hypot(met_gev(chunk["core_mpx"]), met_gev(chunk["core_mpy"])))
            else:
                for label, field in zip(labels, fields, strict=True):
                    fill(histogram, label, met_gev(chunk[field]))
        entries = original["analysis"].num_entries
    save(histogram, labels, entries, output_file)


if __name__ == "__main__":
    app()
