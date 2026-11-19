"""Check that original PHYSLITE MET association data can be read with uproot."""
from __future__ import annotations

import json
from pathlib import Path
from struct import Struct

import awkward as ak
import numpy as np
import typer
import uproot

app = typer.Typer(add_completion=False)

ASSOCIATION = "METAssoc_AnalysisMETAux."
MEMBERWISE_HEADER = Struct(">IHHII")
KEY_AND_INDEX_BYTES = 8
LINK_FIELDS = ("m_persKey", "m_persIndex")


def association_jet_links(branch, n_events: int) -> list:
    """Read split-vector memberwise storage (keys block, then indices block)."""
    decoded = []
    for index in range(branch.num_baskets):
        basket = branch.basket(index)
        for start, stop in zip(basket.byte_offsets, basket.byte_offsets[1:]):
            row = memoryview(basket.data[start:stop])
            encoded_bytes, version, _, _, count = MEMBERWISE_HEADER.unpack_from(row)
            data_size = MEMBERWISE_HEADER.size + count * KEY_AND_INDEX_BYTES
            if (
                not version & 0x4000
                or len(row) != data_size
                or encoded_bytes != (0x40000000 | (data_size - 4))
            ):
                raise ValueError(f"Unrecognized memberwise jetLink basket row {len(decoded)}")
            keys, indices = np.frombuffer(
                row, dtype=">u4", offset=MEMBERWISE_HEADER.size
            ).reshape(2, count)
            decoded.append(list(zip(map(int, keys), map(int, indices))))
            if len(decoded) == n_events:
                return decoded
    raise ValueError(f"Read {len(decoded)} event rows instead of {n_events}")


def raw_link(record: dict) -> tuple[int, int]:
    return record[LINK_FIELDS[0]], record[LINK_FIELDS[1]]


def expected_link(record: dict) -> tuple[int, int]:
    return record["persistent_key"], record["persistent_index"]


@app.command()
def inspect(input_file: Path, release_observation: Path, output_file: Path) -> None:
    oracle = json.loads(release_observation.read_text())
    n_events = len(oracle["events"])
    report = {"input": {}, "association_leaves": {}, "absent_fields": [], "links": {}}
    with uproot.open(input_file, num_workers=1, num_fallback_workers=1) as source:
        tree = source["CollectionTree"]
        report["input"] = {
            "uuid": str(source.file.uuid),
            "events": tree.num_entries,
            "bytes": input_file.stat().st_size,
            "uproot": uproot.__version__,
            "awkward": ak.__version__,
        }
        assert report["input"]["uuid"] == oracle["input"]["uuid"]
        assert report["input"]["events"] == oracle["input"]["entries"]
        assert report["input"]["bytes"] == oracle["input"]["bytes"]
        association = tree[ASSOCIATION]
        leaves = {}
        for name, branch in association.items(recursive=False):
            field = name.rsplit(".", 1)[-1]
            if "xAOD::AuxContainerBase" in name:
                continue
            if field == "jetLink":
                leaves[field] = association_jet_links(branch, n_events)
                shape = "memberwise blocks of uint32 key and index"
            else:
                array = branch.array(entry_stop=n_events, library="ak")
                leaves[field] = ak.to_list(array)
                shape = str(ak.type(array))
            report["association_leaves"][field] = {
                "typename": branch.typename,
                "interpretation": str(branch.interpretation),
                "observed_shape": shape,
            }

        flags = leaves["isMisc"]
        keys_and_indices = leaves["jetLink"]
        objects = leaves["objectLinks"]
        default_jet_links = ak.to_list(
            association["METAssoc_AnalysisMETAux.jetLink"].array(entry_stop=n_events)
        )
        counts = {"jet_links": 0, "object_links": 0, "invalid_object_links": 0}
        name_map = {}
        failing_default_events = 0
        for event_index, event in enumerate(oracle["events"]):
            expected = event["associations"]
            assert event_index == event["entry"]
            assert all(len(rows[event_index]) == len(expected) for rows in leaves.values())
            assert [bool(x) for x in flags[event_index]] == [a["is_misc"] for a in expected]
            expected_jets = [expected_link(row["jet"]) for row in expected]
            assert keys_and_indices[event_index] == expected_jets
            failing_default_events += list(map(raw_link, default_jet_links[event_index])) != expected_jets
            counts["jet_links"] += len(expected_jets)
            for index, row in enumerate(expected):
                expected_objects = row["objects"]
                assert list(map(raw_link, objects[event_index][index])) == [
                    expected_link(link) for link in expected_objects
                ]
                counts["object_links"] += len(expected_objects)
                counts["invalid_object_links"] += sum(not link["valid"] for link in expected_objects)
                for link in [row["jet"], *expected_objects]:
                    if link["valid"]:
                        key = link["persistent_key"]
                        name = link["container"]
                        assert name_map.setdefault(key, name) == name
        report["links"] = {
            **counts,
            "source_container_names_from_root": name_map,
            "events_compared_to_root": n_events,
            "default_jet_link_rows_mismatching_root": failing_default_events,
            "memberwise_jet_link_rows_mismatching_root": 0,
            "object_link_rows_mismatching_root": 0,
        }

        for field in [
            "AnalysisJetsAuxDyn.constituentLinks",
            "AnalysisJetsAuxDyn.originalObjectLink",
            "AnalysisJetsAuxDyn.NNJvt",
            "MET_Core_AnalysisMETAuxDyn.softConstituents",
        ]:
            if field not in tree:
                report["absent_fields"].append(field)
        report["core_and_truth_fields"] = {}
        for stem, expected_key in [
            ("MET_Core_AnalysisMETAuxDyn", "core_terms"),
            ("MET_TruthAuxDyn", "truth_terms"),
        ]:
            fields = {}
            for field in ("mpx", "mpy", "sumet", "source", "name"):
                branch = tree[f"{stem}.{field}"]
                array = branch.array(entry_stop=n_events, library="ak")
                expected_names = [event[expected_key] for event in oracle["events"]]
                if field == "name":
                    assert ak.to_list(array) == expected_names
                else:
                    assert ak.to_list(ak.num(array)) == [len(row) for row in expected_names]
                fields[field] = str(ak.type(array))
            report["core_and_truth_fields"][stem] = fields

    with output_file.open("x") as output:
        json.dump(report, output, indent=2, allow_nan=False)
        output.write("\n")


if __name__ == "__main__":
    app()
