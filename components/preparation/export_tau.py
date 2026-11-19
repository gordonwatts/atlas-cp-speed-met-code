#!/usr/bin/env python3
"""Compile the frozen nominal TauSmearingTool and Tight calibration inputs."""

import hashlib
import json
import math
import os
from pathlib import Path
from typing import Annotated

import typer

RECIPE = "AnalysisBase-25.2.117-mc23-TauTight-NOSYS"
SCALE = "TauAnalysisTools/00-04-00/Smearing/TES_TrueHadTau_RNN_mc23_v3.root"
CALO = "tauRecTools/R22_preprod/CombinedTES_R22_Round2.5_v2.root"
SELECTION = "TauAnalysisAlgorithms/tau_selection_tight_noeleid.conf"
SMOOTHING = "(x<50.)*1.0+((x>=50.)&&(x<70))*(1-(x-50)/20.)"
MODES = ("1p0n", "1p1n", "1pXn", "3p0n", "3pXn")
REALS = ("pt", "eta", "phi", "mass", "charge", "pt_final", "eta_final", "pt_tes", "eta_tes")
INTS = ("charged_tracks", "is_tau_flags", "pan_tau_mode", "has_tes_pt", "tes_compatibility",
        "truth_pdg_id", "truth_is_hadronic", "truth_jet_link_available")


def check(props, expected):
    for key, value in expected.items():
        if str(props.get(key)).strip("'") != str(value):
            raise ValueError(f"Unsupported tau property {key}={props.get(key)!r}, expected {value!r}")


def validate_configuration(path):
    config = json.loads(path.read_bytes())
    smearing = config["seq1_TauSmearingAlg_AnaTauJets"]
    selection = config["seq1_TauSelectionAlg_AnaTauJets_tight"]
    if smearing["Type"] != "CP::TauSmearingAlg" or selection["Type"] != "CP::AsgSelectionAlg":
        raise ValueError("Expected the stock tau preparation algorithms")
    check(smearing["Properties"], {"preselection": ""})
    check(selection["Properties"], {"preselection": "", "selectionDecoration": "selected_tau_tight,as_char"})
    tool = smearing["Properties"]["Tools"]["smearingTool"]
    if tool["Type"] != "TauAnalysisTools::TauSmearingTool":
        raise ValueError("Expected TauSmearingTool")
    check(tool["Properties"], {
        "ApplyFading": 1, "ApplyInsituCorrection": 1, "Campaign": "mc23",
        "InputFilePath": SCALE, "MVATESQualityCheck": 1, "RecommendationTag": "2025-prerec",
        "SkipTruthMatchCheck": 0, "Tool": "", "useFastSim": 0, "useGNTau": 0,
    })
    internal = tool["Properties"]["Tools"]["TauSmearingTool"]["Properties"]
    check(internal, {"ApplyFading": 1, "ApplyInsituCorrection": 1, "InputFilePath": SCALE,
                     "MVATESQualityCheck": 1, "SkipTruthMatchCheck": 0})
    check(internal["Tools"]["TauResolutionTool"]["Properties"], {"WeightFileName": CALO})
    cut_tool = selection["Properties"]["Tools"]["selectionTool"]
    if cut_tool["Type"] != "TauAnalysisTools::TauSelectionTool":
        raise ValueError("Expected TauSelectionTool")
    check(cut_tool["Properties"], {
        "ConfigPath": SELECTION, "AbsCharge": 1, "PtMin": 20, "SelectionCuts": 187,
        "EleIDWP": 1, "EleIDVersion": 1, "JetIDWP": 9, "CreateControlPlots": 0, "MuonOLR": 0,
    })


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def numbers(values):
    return " ".join(repr(float(value)) for value in values)


def axis_line(axis):
    count = axis.GetNbins()
    bins = axis.GetXbins()
    uniform = bins.GetSize() == 0
    edges = [float(axis.GetBinLowEdge(i)) for i in range(1, count + 2)] if uniform else [
        float(bins.At(i)) for i in range(bins.GetSize())
    ]
    if uniform:
        edges[0], edges[-1] = float(axis.GetXmin()), float(axis.GetXmax())
    if len(edges) != count + 1 or any(not math.isfinite(x) for x in edges):
        raise ValueError("Unknown tau scale axis")
    return f"{int(uniform)} {len(edges)} {numbers(edges)}"


def export_scale(root, path):
    source = root.TFile.Open(str(path), "READ")
    if not source or source.IsZombie():
        raise RuntimeError(f"Invalid tau scale factor file {path}")
    try:
        calibrations = {}
        locations = {}
        absolute_eta = False
        # Match ReadInputs ordering, including repeated sf keys.
        for key in source.GetListOfKeys():
            name = str(key.GetName())
            if name == "Yaxis" and str(key.ReadObj().GetTitle()) == "|eta|":
                absolute_eta = True
            if name not in ("sf_1p", "sf_3p"):
                continue
            hist = key.ReadObj()
            nx, ny = hist.GetNbinsX(), hist.GetNbinsY()
            if str(hist.ClassName()) != "TH2F" or hist.GetNcells() != (nx + 2) * (ny + 2):
                raise ValueError(f"Unsupported histogram shape {name}")
            title = str(hist.GetTitle())
            if "".join(title.split()) != SMOOTHING:
                raise ValueError(f"Unexpected nominal fade for {name}: {title}")
            calibrations[name] = "\n".join([
                "TauNominalScaleV1", axis_line(hist.GetXaxis()), axis_line(hist.GetYaxis()),
                numbers(hist.GetBinContent(i) for i in range(hist.GetNcells())),
            ])
            locations[name] = {"cycle": key.GetCycle(), "class": str(hist.ClassName()),
                               "bins": [nx, ny], "title": title}
        if not absolute_eta or calibrations.keys() != {"sf_1p", "sf_3p"}:
            raise ValueError("Missing absolute eta tau scale calibration")
        return calibrations, locations
    finally:
        source.Close()


def export_calo(root, path):
    source = root.TFile.Open(str(path), "READ")
    if not source or source.IsZombie():
        raise RuntimeError(f"Invalid combined TES file {path}")
    try:
        rows = ["METSignificanceTauV1"]  # Same graph record contract as the significance reader.
        details = []
        for mode in MODES:
            for bin_index in range(5):
                name = f"tauRec/Graph_from_ResolutionEt_tauRec_{mode}_{bin_index}"
                graph = source.Get(name)
                if not graph or not graph.InheritsFrom("TGraph"):
                    raise ValueError(f"Missing calo-only resolution graph {name}")
                size = graph.GetN()
                abscissa = [float(graph.GetPointX(i)) for i in range(size)]
                ordinate = [float(graph.GetPointY(i)) for i in range(size)]
                if not abscissa or any(not math.isfinite(x) for x in abscissa):
                    raise ValueError(f"Invalid tau calo graph abscissa {name}")
                sorted_x = bool(graph.TestBit(root.TGraph.kIsSortedX))
                maximum = max(abscissa)
                rows += [f"graph {int(sorted_x)} {maximum!r} {size}", numbers(abscissa), numbers(ordinate)]
                details.append({"path": name, "class": str(graph.ClassName()), "points": size,
                                "sorted_x": sorted_x, "max_pt_mev": maximum})
        return "\n".join(rows) + "\n", details
    finally:
        source.Close()


def export_cuts(root, path):
    settings = root.TEnv(str(path))
    cuts = {key: str(settings.GetValue(key, "")) for key in
            ("SelectionCuts", "PtMin", "AbsEtaRegion", "AbsCharge", "NTracks", "JetIDWP", "EleIDWP")}
    if cuts["SelectionCuts"].split() != ["PtMin", "AbsEtaRegion", "AbsCharge", "NTracks", "JetIDWP", "EleIDWP"]:
        raise ValueError(f"Unexpected tau selection flags {cuts['SelectionCuts']}")
    check(cuts, {"PtMin": 20, "AbsCharge": 1, "JetIDWP": "JETIDRNNTIGHT", "EleIDWP": "ELEIDNONE"})
    if [int(x) for x in cuts["NTracks"].split(";")] != [1, 3]:
        raise ValueError(f"Unexpected tau track multiplicities {cuts['NTracks']}")
    edges = [float(x) for x in cuts["AbsEtaRegion"].split(";")]
    if edges != [0., 1.37, 1.52, 2.5]:
        raise ValueError(f"Unexpected tau eta regions {edges}")
    return cuts, dict(zip(("eta_low", "eta_crack_low", "eta_crack_high", "eta_high"), edges))


def extension(name, inputs, config, shape=None):
    node = {"nodetype": "extension", "type": f"atlas.met.{name}", "version": 1,
            "children": [], "config": config}
    output = {"name": "values", "type": "real"}
    if shape:
        node["output_shape"] = [len(shape)]
        output.update(shape=[len(shape)], description=", ".join(shape))
    return {"name": name, "version": 1, "data": node, "output": output,
            "inputs": [{"name": name, "type": kind} for name, kind in inputs]}


def write_json(path, obj):
    with path.open("x", encoding="utf-8") as destination:
        json.dump(obj, destination, separators=(",", ":"), allow_nan=False)
        destination.write("\n")


def main(
    reference_config: Annotated[Path, typer.Option(exists=True, dir_okay=False)],
    output_dir: Annotated[Path, typer.Option()],
):
    """Generate tau.json and its source receipt in a fresh directory."""
    import ROOT

    if os.environ.get("AtlasVersion") != "25.2.117" or str(ROOT.gROOT.GetVersion()).replace("/", ".") != "6.40.02":
        raise RuntimeError("Source AnalysisBase 25.2.117 first")
    if ROOT.gSystem.Load("libPathResolverDict") < 0:
        raise RuntimeError("Missing release PathResolver")
    validate_configuration(reference_config)
    if output_dir.exists():
        raise FileExistsError(output_dir)
    paths = {}
    for logical in (SCALE, CALO, SELECTION):
        resolved = ROOT.PathResolver.FindCalibFile(logical)
        if not resolved:
            raise RuntimeError(f"Could not resolve {logical}")
        paths[logical] = Path(str(resolved)).resolve(strict=True)
    hist, details = export_scale(ROOT, paths[SCALE])
    graph, graphs = export_calo(ROOT, paths[CALO])
    cuts, eta = export_cuts(ROOT, paths[SELECTION])
    config = {"recipe": RECIPE, "calo_graphs": graph, "min_pt_gev": float(cuts["PtMin"]), **hist, **eta}
    shapes = ("pt", "eta", "phi", "mass", "selected_tau_tight", "reached_cut_pass_bits",
              "TESCompatibility", "smear_branch", "scale")
    payload = {"schema_version": 3, "corrections": [
        extension("tau_track_charged", [("flag_set", "int")], {}),
        extension("tau_prepare", [(name, "real") for name in REALS] + [(name, "int") for name in INTS], config, shapes),
    ]}
    output_dir.mkdir(parents=True, exist_ok=False)
    write_json(output_dir / "tau.json", payload)
    write_json(output_dir / "provenance.json", {
        "athena": "6ce3f7ae7c41dd9eaec170c3c349df15faef7e5c",
        "analysisbase": "25.2.117", "root": str(ROOT.gROOT.GetVersion()), "plugin_api": 4,
        "correctionlib": "3dcfa20ec01fb40c4d5fec58383785519dec1c00", "recipe": RECIPE,
        "reference_config": str(reference_config.resolve()), "config_sha256": digest(reference_config),
        "sources": {logical: {"resolved": str(path), "sha256": digest(path)} for logical, path in paths.items()},
        "nominal_histograms": details, "calo_graphs": graphs, "cuts": cuts, "outputs": shapes,
        "tau_json_sha256": digest(output_dir / "tau.json"),
    })


if __name__ == "__main__":
    typer.run(main)
