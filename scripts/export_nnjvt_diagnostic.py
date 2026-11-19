"""Wrap the already implemented and independently identified upstream NNJVT nodes."""
from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path

import typer

app = typer.Typer(add_completion=False)


def load_exporter(path: Path):
    spec = importlib.util.spec_from_file_location("pinned_small_jvt_exporter", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@app.command()
def create(
    prior_exporter: Path, group_data: Path, full_config: Path, output: Path,
) -> None:
    exporter = load_exporter(prior_exporter)
    config = json.loads(full_config.read_text())
    maker = config["seq1_MetMakerAlg_AnaMET"]["Properties"]["Tools"]["makerTool"]["Properties"]
    jvt = maker["Tools"]["JvtSelTool"]["Properties"]
    expected_config = {
        "ConfigDir": "JetPileupTag/NNJvt/2022-03-22",
        "ConfigFile": "NNJVT.Cuts.FixedEffPt.Offline.Nonprompt_All_MaxW.json",
        "WorkingPoint": "FixedEffPt",
        "JetEtaName": "eta",
    }
    for name, value in expected_config.items():
        assert jvt[name].strip("'") == value
    assert float(jvt["MinPtForJvt"]) == exporter.WINDOW["min_pt"]
    assert float(jvt["MaxPtForJvt"]) == exporter.WINDOW["max_pt"]
    assert float(jvt["MinEtaForJvt"]) == exporter.WINDOW["min_abs_eta"]
    assert float(jvt["MaxEtaForJvt"]) == exporter.WINDOW["max_abs_eta"]
    files = {}
    for name in ("network", "cuts"):
        relative, expected = exporter.SOURCE_FILES[name]
        path = group_data / relative
        raw = path.read_bytes()
        assert hashlib.sha256(raw).hexdigest() == expected, path
        files[name] = raw
    cut_config = {"cut_map": exporter.encode_cut_map(files["cuts"]), "mev_per_gev": 1000.0}
    pt_eta = [
        exporter.variable("pt", "real", "Float-stored pT in MeV."),
        exporter.variable("eta", "real", "Float-stored signed eta."),
    ]
    update = exporter.correction(
        "small_nnjvt_update", "nnjvt_update",
        pt_eta + [
            exporter.variable("sum_pt_trk", "real", "Primary-vertex track sum in MeV."),
            exporter.variable("jvf_corr", "real", "Stored JVFCorr."),
            exporter.variable("vertex_kind", "int", "Release primary-vertex classification."),
            exporter.variable("feature_status", "int", "Availability of input moments."),
        ],
        {**cut_config, "max_pt": exporter.WINDOW["max_pt"], "network_json": files["network"].decode("utf-8")},
        "Nominal NNJVT update for MET preparation.",
    )
    selection = exporter.correction(
        "small_nnjvt_select", "nnjvt_select",
        pt_eta + [
            exporter.variable("score", "real", "New float NNJVT score."),
            exporter.variable("score_present", "int", "Availability of the new score."),
        ],
        {**exporter.WINDOW, **cut_config},
        "MET FixedEffPt jet selector.",
    )
    with output.open("x") as file:
        json.dump({"schema_version": 3, "description": "MET NNJVT diagnostic", "corrections": [update, selection]}, file)
        file.write("\n")


if __name__ == "__main__":
    app()
