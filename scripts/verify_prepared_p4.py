"""Match read-only prepared-particle getter records to the untouched stock output."""
from __future__ import annotations

import json
from pathlib import Path

import awkward as ak
import numpy as np
import typer

from histogram_reference import root_file
from summarize_observer import COLLECTIONS, MAKER_LEPTONS
from verify_status_trace import check_stock_parity

app = typer.Typer(add_completion=False)


@app.command()
def check(input_file: Path, original: Path, traced: Path, output: Path) -> None:
    rows = {prefix: [] for prefix in COLLECTIONS}
    with (traced / "prepared-p4.tsv").open() as file:
        for line in file:
            entry, prefix, row, key, index, passed, *four_vector = line.rstrip("\n").split("\t")
            rows[prefix].append(
                [int(entry), int(row), int(key), int(index), int(passed), *map(float.fromhex, four_vector)]
            )
    config = json.loads((original / "full_config.json").read_text())
    maker = config["seq1_MetMakerAlg_AnaMET"]["Properties"]
    selection_expression = {
        prefix: maker[f"{species}Selection"] for species, prefix in MAKER_LEPTONS.items()
    }
    with (
        root_file(input_file) as raw,
        root_file(original / "reference.root") as original_file,
        root_file(traced / "reference.root") as traced_file,
    ):
        a = original_file["analysis"]
        candidate = traced_file["analysis"]
        events = candidate.num_entries
        branch_count = check_stock_parity(a, candidate)
        report = {"events": events, "stock_branches_unchanged": branch_count, "collections": {}}
        for prefix, collection in COLLECTIONS.items():
            parts = np.array(rows[prefix], dtype=np.float64)
            key = candidate[f"{prefix}_originalObjectLink"][f"{prefix}_originalObjectLink.m_persKey"].array()
            index = candidate[f"{prefix}_originalObjectLink"][f"{prefix}_originalObjectLink.m_persIndex"].array()
            counts = ak.to_numpy(ak.num(index))
            np.testing.assert_array_equal(np.repeat(np.arange(events), counts), parts[:, 0])
            np.testing.assert_array_equal(ak.to_numpy(ak.flatten(ak.local_index(index))), parts[:, 1])
            np.testing.assert_array_equal(ak.to_numpy(ak.flatten(key)), parts[:, 2])
            np.testing.assert_array_equal(ak.to_numpy(ak.flatten(index)), parts[:, 3])
            counts_wp = None
            selected = None
            if prefix in selection_expression:
                gate = ak.ones_like(index, dtype=bool)
                for item in selection_expression[prefix].split("&&"):
                    aux, _, mode = item.partition(",")
                    step = candidate[f"{prefix}_{aux}"].array()
                    gate = gate & ((step != 0) if mode == "as_char" else (step == 0xffffffff))
                selected = ak.to_numpy(ak.flatten(gate))
                np.testing.assert_array_equal(selected, parts[:, 4])
                counts_wp = int(np.count_nonzero(selected))
            else:
                assert np.all(parts[:, 4] == -1)
            has_raw_mass = f"{collection}AuxDyn.m" in raw["CollectionTree"]
            for axis, field in enumerate(("pt", "eta", "phi", "m"), start=5):
                if field == "m" and not has_raw_mass:
                    continue
                expected = ak.to_numpy(ak.flatten(candidate[f"{prefix}_{field}"].array()))
                np.testing.assert_array_equal(parts[:, axis].astype(expected.dtype), expected)
            mass = parts[:, 8]
            energy = parts[:, 9]
            report["collections"][prefix] = {
                "all_rows": len(parts), "selected_rows": counts_wp, "original_mass_leaf_exists": has_raw_mass,
                "actual_mass_getter_MeV_min_max": [float(mass.min()), float(mass.max())],
                "actual_energy_getter_MeV_min_max": [float(energy.min()), float(energy.max())],
                "selected_mass_getter_MeV_min_max": [float(mass[selected].min()), float(mass[selected].max())] if counts_wp else None,
                "selected_energy_getter_MeV_min_max": [float(energy[selected].min()), float(energy[selected].max())] if counts_wp else None,
                "actual_nonfinite_getter_values": int(np.count_nonzero(~np.isfinite(parts[:, 5:]))),
            }
    with output.open("x") as file:
        json.dump(report, file, indent=2, allow_nan=False)
        file.write("\n")


if __name__ == "__main__":
    app()
