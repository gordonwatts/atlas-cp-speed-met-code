"""Release-independent control reading persisted PVSoftTrkCore from PHYSLITE."""
from pathlib import Path

import awkward as ak
import numpy as np
import typer

from histogram_reference import book, fill, met_gev, root_file, save

app = typer.Typer(add_completion=False)
LABEL = "PVSoftTrkCore"
STEM = "MET_Core_AnalysisMETAuxDyn."


@app.command()
def read(
    input_file: Path,
    output_file: Path,
    events: int = typer.Option(20000, min=1),
    bin_output: Path | None = None,
) -> None:
    histogram = book([LABEL])
    fields = [STEM + f for f in ("name", "mpx", "mpy")]
    entries = 0
    ordered_bins = []
    with root_file(input_file) as original:
        tree = original["CollectionTree"]
        for chunk in tree.iterate(fields, entry_stop=events, step_size=2000, library="ak"):
            selected = chunk[fields[0]] == LABEL
            if not ak.all(ak.sum(selected, axis=1) == 1):
                raise ValueError("Expected one persisted PVSoftTrkCore in every entry.")
            x = ak.to_numpy(ak.firsts(chunk[fields[1]][selected], axis=1))
            y = ak.to_numpy(ak.firsts(chunk[fields[2]][selected], axis=1))
            gev = np.hypot(met_gev(x), met_gev(y))
            fill(histogram, LABEL, gev)
            if bin_output is not None:
                ordered_bins.append(histogram.axes[1].index(gev) + 1)
            entries += len(chunk)
    if entries != events:
        raise ValueError(f"Read {entries} entries, expected {events}.")
    save(histogram, [LABEL], entries, output_file)
    if bin_output is not None:
        with bin_output.open("xb") as output:
            np.savez(output, entry=np.arange(entries), flow_inclusive_bin=np.concatenate(ordered_bins))


if __name__ == "__main__":
    app()
