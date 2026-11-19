"""Record measured fixed-input invariants and source exceptions before reuse trials."""
from pathlib import Path
import json

import typer

from summarize_observer import FINAL_QUANTITIES, TERMS

app = typer.Typer(add_completion=False)
HARD_QUANTITIES = ("mpx", "mpy", "sumet", "met", "phi")


@app.command()
def record(stock: Path, stock_status: Path, candidate: Path, output: Path) -> None:
    observed = json.loads((stock / "summary.json").read_text())
    traced = json.loads((stock_status / "verification.json").read_text())
    check = json.loads((candidate / "comparison.json").read_text())
    decisions = json.loads((candidate / "trace-comparison.json").read_text())
    n_events = observed["entries"]
    assert n_events == traced["entries"] == check["events"] == decisions["entries"] == 20000
    assert traced["trace_rows_aligned_with_entries_and_states"] and check["passed"] and decisions["passed"]
    assert not check["discrete_failures"]
    config = json.loads((stock / "full_config.json").read_text())
    properties = config["seq1_MetSignificanceAlg_AnaMET"]["Properties"]["Tools"]["significanceTool"]["Properties"]
    hard = [term for term in TERMS if term not in ("PVSoftTrk", "Final")]
    report = {"events": n_events, "states": {}}
    for label, state in observed["states"].items():
        if label != "NOSYS":
            assert all(state["changed_from_nominal"][field][term] == 0 for field in HARD_QUANTITIES for term in hard)
        counters = traced["states"][label]["calls_by_action"]
        expected_draws = 0 if label == "NOSYS" else n_events
        assert counters.get("gaussian", 0) == expected_draws
        report["states"][label] = {
            "stock_hard_terms_exact_vs_nominal": hard,
            "stock_ordered_constituent_links_and_weights_exact_vs_nominal": True,
            "stock_association_and_jet_jvt_sequences_exact_vs_nominal": True,
            "stock_seed_exact_vs_nominal": True,
            "stock_first_gaussian_exact_across_active_states": label != "NOSYS",
            "stock_gaussian_calls": expected_draws,
            "stock_significance_changes_vs_nominal": {
                quantity: state["changed_from_nominal"].get(quantity, {}).get("Final", 0)
                for quantity in FINAL_QUANTITIES
            },
            "association_writes": counters["associationFlag"],
            "jvt_decisions": counters["jetJvt"],
        }
    report["B_ordered_stock_decision_agreement"] = decisions["exact_ordered_counts"]
    report["B_scalar_outputs_exact_to_A"] = all(row["exact_finite"] == row["values"] for row in check["numerical_checks"].values())
    report["covariance"] = {
        "observed_stock_soft_resolution": float(properties["SoftTermReso"]),
        "stock_soft_term_model_parameter": int(properties["SoftTermParam"]),
        "observed_complete_valid_final_and_soft_rows": n_events * len(observed["states"]),
        "dependence": "Longitudinal, transverse and mixed covariance projections use the Final MET direction in each state.",
        "state_carry_source_conditions": [
            "Wrong total-source mask returns before recomputing the significance.",
            "When longitudinal variance equals zero, the release retains its prior squared significance and rho.",
            "Within the matched valid sample the fixed track-soft term contributes nonzero isotropic variance; retain the release state across calls for the distinct sentinel paths.",
        ],
        "source": "AnalysisBase 25.2.117 METUtilities/Root/METSignificance.cxx varianceMET and AddSoftTerm",
    }
    report["raw_authority"] = {
        "A": str(stock.resolve()), "A_status": str(stock_status.resolve()), "B": str(candidate.resolve()),
    }
    with output.open("x") as file:
        json.dump(report, file, indent=2, allow_nan=False)
        file.write("\n")


if __name__ == "__main__":
    app()
