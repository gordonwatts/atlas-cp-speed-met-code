"""Validate a stock MET observation and summarize discrete/numerical coverage."""
from __future__ import annotations

import json
from pathlib import Path

import awkward as ak
import numpy as np
import typer
import uproot

app = typer.Typer(add_completion=False)
COLLECTIONS = {
    "jet": "AnalysisJets", "ele": "AnalysisElectrons", "pho": "AnalysisPhotons",
    "mu": "AnalysisMuons", "tau": "AnalysisTauJets",
}
TERMS = {
    "RefEle": ("ele", 65537), "RefGamma": ("pho", 65538), "RefTau": ("tau", 65540),
    "Muons": ("mu", 65544), "MuonEloss": ("mu", 262152), "RefJet": ("jet", 65552),
    "PVSoftTrk": ("jet", 81952), "Final": (None, 65600),
}
MAKER_LEPTONS = {"electrons": "ele", "photons": "pho", "muons": "mu", "taus": "tau"}
FINAL_QUANTITIES = ("significance", "sigDirectional", "METOverSqrtSumET", "METOverSqrtHT")
TERM_QUANTITIES = ("mpx", "mpy", "sumet", "met", "phi") + FINAL_QUANTITIES


def difference(before, after, mask) -> dict:
    different = before != after
    if mask is not None:
        different_selected = different[mask]
        selected_changed = int(ak.sum(different_selected))
    else:
        selected_changed = None
    return {
        "changed": int(ak.sum(different)),
        "changed_selected": selected_changed,
        "max_absolute_difference": float(ak.max(abs(before - after), initial=0)),
    }


def upstream(source, observed, maker, n_events: int) -> tuple:
    selection_expression = {
        prefix: maker[f"{species}Selection"] for species, prefix in MAKER_LEPTONS.items()
    }
    report = {}
    keys_by_prefix = {}
    original_indices = {}
    masks = {}
    for prefix, collection in COLLECTIONS.items():
        link_branch = observed[f"{prefix}_originalObjectLink"]
        key = link_branch[f"{prefix}_originalObjectLink.m_persKey"].array()
        index = link_branch[f"{prefix}_originalObjectLink.m_persIndex"].array()
        raw_base = collection + "AuxDyn."
        raw_pt = source[raw_base + "pt"].array(entry_stop=n_events)
        assert ak.all(ak.num(index) == ak.num(raw_pt))
        assert ak.all(index == ak.local_index(raw_pt))
        unique_keys = np.unique(ak.to_numpy(ak.flatten(key))).tolist()
        assert len(unique_keys) <= 1
        keys_by_prefix[prefix] = unique_keys
        original_indices[prefix] = ak.to_list(index)
        row = {"object_count": int(ak.sum(ak.num(index))), "original_link_keys": unique_keys}
        if prefix in selection_expression:
            row["selection_expression"] = selection_expression[prefix]
            mask = ak.ones_like(index, dtype=bool)
            row["gate_passing"] = {}
            for token in selection_expression[prefix].split("&&"):
                aux, _, accessor = token.partition(",")
                value = observed[f"{prefix}_{aux}"].array()
                assert ak.all(ak.num(value) == ak.num(index))
                gate = (value != 0) if accessor == "as_char" else (value == 0xffffffff)
                row["gate_passing"][token] = int(ak.sum(gate))
                mask = mask & gate
            row["passing_working_point"] = int(ak.sum(mask))
            masks[prefix] = ak.to_list(mask)
        else:
            mask = None
            old = source[raw_base + "NNJvtPass"].array(entry_stop=n_events)
            new = observed["jet_NNJvtPass"].array()
            row["recomputed_NNJvtPass_mismatches"] = int(ak.sum(old != new))
            row["recomputed_NNJvtPass_counts"] = int(ak.sum(new != 0))
            row["stored_NNJvt_score_available"] = raw_base + "NNJvt" in source
        row["original_vs_nominal"] = {}
        for field in ("pt", "eta", "phi", "m"):
            if raw_base + field not in source:
                row["original_vs_nominal"][field] = {"raw_aux_present": False}
                continue
            before = source[raw_base + field].array(entry_stop=n_events)
            after = observed[f"{prefix}_{field}"].array()
            row["original_vs_nominal"][field] = difference(before, after, mask)
        report[prefix] = row
    return report, keys_by_prefix, original_indices, masks


def check_links(links, weights, original_indices, masks, prepared_keys) -> dict:
    counts = dict.fromkeys(TERMS, 0)
    hard = set(TERMS) - {"PVSoftTrk", "Final"}
    for event, terms in enumerate(ak.to_list(links)):
        weights_by_term = ak.to_list(weights[event])
        for index, (term, (prefix, _)) in enumerate(TERMS.items()):
            row = terms[index]
            if term in hard:
                assert len(row) == len(weights_by_term[index])
            elif term == "Final":
                assert not row and not weights_by_term[index]
            counts[term] += len(row)
            for link in row:
                position = link["m_persIndex"]
                key = link["m_persKey"]
                assert position < len(original_indices[prefix][event])
                assert prepared_keys.setdefault(prefix, key) == key
                if term in ("RefEle", "RefGamma", "RefTau", "Muons"):
                    assert masks[prefix][event][position]
    return counts


def check_states(observed, states, event_count, original_indices, masks) -> dict:
    baseline = {}
    report = {}
    prepared_keys = {}
    for state in states:
        arrays = {
            field: observed[f"fullMet_{field}_{state}"].array()
            for field in TERM_QUANTITIES + ("term", "source", "constituentLinks", "constituentWeights")
        }
        assert all(len(rows) == event_count for rows in arrays.values())
        assert ak.to_list(arrays["term"]) == [list(TERMS)] * event_count
        assert ak.to_list(arrays["source"]) == [[mask for _, mask in TERMS.values()]] * event_count
        report[state] = {"invalid_or_nonfinite": {}, "changed_from_nominal": {}}
        for field in TERM_QUANTITIES:
            data = arrays[field]
            assert ak.all(ak.num(data) == len(TERMS))
            report[state]["invalid_or_nonfinite"][field] = int(ak.sum(~np.isfinite(data)))
            if field in FINAL_QUANTITIES or field in ("met", "phi", "sumet"):
                final = observed[f"finalMet_{field}_{state}"].array()
                np.testing.assert_array_equal(ak.to_numpy(final), ak.to_numpy(data[:, -1]))
            if state == "NOSYS":
                continue
            changes = data != baseline[field]
            report[state]["changed_from_nominal"][field] = {
                term: int(ak.sum(changes[:, index])) for index, term in enumerate(TERMS)
            }
            if field in ("mpx", "mpy", "sumet", "met", "phi"):
                np.testing.assert_array_equal(ak.to_numpy(data[:, :-2]), ak.to_numpy(baseline[field][:, :-2]))
        if state == "NOSYS":
            counts = check_links(
                arrays["constituentLinks"], arrays["constituentWeights"], original_indices, masks, prepared_keys
            )
            report[state]["counts_by_term"] = counts
            baseline = arrays
        else:
            for field in ("constituentLinks", "constituentWeights"):
                assert ak.to_list(arrays[field]) == ak.to_list(baseline[field])
            assert any(report[state]["changed_from_nominal"]["met"][name] for name in ("PVSoftTrk", "Final"))
    return {"states": report, "prepared_link_keys": prepared_keys}


def compare_nominal(left, right, fields) -> int:
    for field in fields:
        left_values = left[field].array()
        right_values = right[field].array(entry_stop=len(left_values))
        try:
            np.testing.assert_array_equal(ak.to_numpy(ak.flatten(left_values, axis=None)), ak.to_numpy(ak.flatten(right_values, axis=None)))
        except (TypeError, ValueError):
            assert ak.to_list(left_values) == ak.to_list(right_values), field
    return len(fields)


@app.command()
def summarize(
    input_file: Path,
    result_path: Path,
    output_file: Path,
    compare: Path | None = None,
) -> None:
    expected = ["NOSYS", "MET_SoftTrk_ResoPara", "MET_SoftTrk_ResoPerp", "MET_SoftTrk_Scale__1down", "MET_SoftTrk_Scale__1up"]
    with uproot.open(result_path / "hist-reference.root", num_workers=1, num_fallback_workers=1) as meta:
        states = meta["systematicsList"].axis().labels()
        assert states in (expected, ["NOSYS"], ["NOSYS", "MET_SoftTrk_ResoCorr"])
        assert "systematicsListOnlyWeights" not in meta
    config = json.loads((result_path / "full_config.json").read_text())
    maker = config["seq1_MetMakerAlg_AnaMET"]["Properties"]
    report = {"result": str(result_path), "state_order": states}
    with (
        uproot.open(input_file, num_workers=1, num_fallback_workers=1) as source_file,
        uproot.open(result_path / "reference.root", num_workers=1, num_fallback_workers=1) as observer_file,
    ):
        source = source_file["CollectionTree"]
        observed = observer_file["analysis"]
        n_events = observed.num_entries
        report["entries"] = n_events
        for field in ("eventNumber", "runNumber", "mcChannelNumber", "averageInteractionsPerCrossing", "actualInteractionsPerCrossing"):
            before = source["EventInfoAuxDyn." + field].array(entry_stop=n_events)
            after = observed[field].array()
            np.testing.assert_array_equal(ak.to_numpy(before), ak.to_numpy(after))
        report["runNumbers"] = np.unique(observed["runNumber"].array(library="np")).tolist()
        report["mcChannelNumbers"] = np.unique(observed["mcChannelNumber"].array(library="np")).tolist()
        report["identity_matches_raw_all_entries"] = True
        upstream_report, _, original_indices, masks = upstream(source, observed, maker, n_events)
        report["upstream"] = upstream_report
        report.update(check_states(observed, states, n_events, original_indices, masks))
        if compare:
            with uproot.open(compare, num_workers=1, num_fallback_workers=1) as baseline_file:
                baseline = baseline_file["analysis"]
                fields = [k for k in baseline.keys() if k in observed and "Link" not in k]
                report["other_observer_common_branches_equal"] = compare_nominal(baseline, observed, fields)
    report["term_order_and_source"] = {term: value[1] for term, value in TERMS.items()}
    report["observed_statuses"] = {
        "cp_job_success": "Completion and matching entry/state/term coverage; ANA_CHECK rejects failed stages.",
        "constituent_absence": "Observed as ordered empty vectors, distinct from skipped/failed function calls.",
        "individual_tool_statuses": "CP MET algorithms consume return codes and do not publish per-call codes.",
        "internal_selection_flags": "The stock ntuple exposes upstream selection decorations, input association links, and output constituent choices. Transient flags inside METMaker require a separate stock observation.",
    }
    with output_file.open("x") as output:
        json.dump(report, output, indent=2, allow_nan=False)
        output.write("\n")


if __name__ == "__main__":
    app()
