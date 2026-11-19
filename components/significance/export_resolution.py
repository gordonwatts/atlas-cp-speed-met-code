#!/usr/bin/env python3
"""Export the fixed MET significance inputs using stock PathResolver and ROOT."""

from contextlib import contextmanager
import hashlib
import json
import math
import os
from pathlib import Path
from typing import Annotated

import typer

ATHENA = "6ce3f7ae7c41dd9eaec170c3c349df15faef7e5c"
CORE = "3dcfa20ec01fb40c4d5fec58383785519dec1c00"
JET_CONFIG = "JES_data2017_2016_2015_Recommendation_PFlow_Aug2018_rel21.config"
JET_AREA = "CalibArea-00-04-81/"
EGAMMA_MODEL = "es2024_Run3_v0"
TAU_FILE = "tauRecTools/R22_preprod/CombinedTES_R22_Round2.5_v2.root"
MUON_RELEASE = "Recs2025_12_11_Run2Run3"
MODES = ("1p0n", "1p1n", "1pXn", "3p0n", "3pXn")


def check_fields(properties, expected):
    for name, value in expected.items():
        actual = str(properties.get(name)).strip("'")
        if actual != str(value):
            raise ValueError(f"Unsupported setting {name}={actual!r}. Expected {value!r}.")


def validate_config(config, algorithm):
    alg = config[algorithm]
    if alg["Type"] != "CP::MetSignificanceAlg":
        raise ValueError("Expected stock MetSignificanceAlg in the reference config")
    settings = alg["Properties"]
    for key, default in (("totalMETName", "Final"), ("jetTermName", "RefJet"),
                         ("softTermName", "PVSoftTrk")):
        if key in settings:
            check_fields(settings, {key: default})
    tool = settings["Tools"]["significanceTool"]
    if tool["Type"] != "met::METSignificance":
        raise ValueError("Expected stock METSignificance tool in the reference config")
    props = tool["Properties"]
    check_fields(props, {
        "ApplyBias": 0, "ScalarBias": 0, "DoJerForEMu": 0, "DoPhiReso": 0,
        "EMuResoAux": "", "EgammaDecorrelationModel": "1NP_v1", "EgammaESModel": EGAMMA_MODEL,
        "EgammaUseFastsim": 0, "IsDataJet": 0, "IsDataMuon": 0, "JetCalibArea": "00-04-81",
        "JetCalibConfig": JET_CONFIG, "JetCalibSequence": "JetArea_Residual_EtaJES_GSC_Smear",
        "JetCollection": "AntiKt4EMPFlow", "JetEtaThr": -1, "JetPtThr": -1, "JetResoAux": "",
        "MuonCalibMode": -1, "SoftTermParam": 0, "SoftTermReso": 8.5, "TauTESConfig": TAU_FILE,
        "TreatPUJets": 0, "egammaCalibTool": "", "jetCalibTool": "", "tauResolutionTool": "",
    })
    children = props["Tools"]
    check_fields(children["METSigAutoConf_EgammaCalibrationAndSmearingTool"]["Properties"], {
        "ESModel": EGAMMA_MODEL, "ResolutionType": "SigmaEff90", "decorrelationModel": "1NP_v1",
        "useFastSim": 0, "useMVACalibration": 0,
    })
    tau = children["METSigAutoConf_TauPerfTool"]["Properties"]
    check_fields(tau, {"WeightFileName": TAU_FILE})
    jet_name = "MetSigAutoConf_JetCalibTool_AntiKt4EMPFlow"
    jet_properties = []
    if jet_name in children:
        jet_properties.append(children[jet_name]["Properties"])
    # In AnalysisBase 25.2.117 standalone the full config writer can flatten
    # this privately created tool into the adjacent tau property block.
    if jet_name in tau:
        prefix, _, name = str(tau[jet_name]).strip("'").partition("/")
        if prefix not in ("18JetCalibrationTool", "JetCalibrationTool") or name != f"{props['Name']}.{jet_name}":
            raise ValueError(f"Unexpected resolved jet tool marker {tau[jet_name]!r}")
        jet_properties.append(tau)
    if not jet_properties:
        raise ValueError("Missing resolved JetCalibrationTool properties in the stock config")
    for jet in jet_properties:
        check_fields(jet, {
            "CalibArea": "00-04-81", "CalibSequence": "JetArea_Residual_EtaJES_GSC_Smear",
            "ConfigDir": "JetCalibTools/CalibrationConfigs/", "ConfigFile": JET_CONFIG,
            "DEVmode": 0, "IsData": 0, "JetCollection": "AntiKt4EMPFlow",
        })
    check_fields(children["MuonCalibTool"]["Properties"], {
        "calibMode": 0, "release": MUON_RELEASE, "systematicScheme": "Corr_Scale",
        "useRandomRunNumber": 0, "expertMode_RunNumber": 0, "expert_validationMode": 0,
        "ScaleAndSmearTool": "",
    })


def sha256(path):
    with path.open("rb") as source:
        return hashlib.file_digest(source, "sha256").hexdigest()


class Sources:
    def __init__(self, root):
        self.root = root
        self.files = {}

    def resolve(self, relative):
        found = self.root.PathResolver.FindCalibFile(relative)
        if not found:
            raise RuntimeError(f"PathResolver could not resolve {relative}")
        path = Path(str(found)).resolve(strict=True)
        self.files[relative] = {
            "path": str(path), "sha256": sha256(path), "size": path.stat().st_size, "objects": {},
        }
        return path

    @contextmanager
    def root_file(self, relative):
        path = self.resolve(relative)
        source = self.root.TFile.Open(str(path), "READ")
        if not source or source.IsZombie():
            raise RuntimeError(f"Cannot open {path}")
        try:
            yield source
        finally:
            source.Close()

    def object(self, source, relative, key, base, dimension=None):
        obj = source.Get(key)
        if not obj or not obj.InheritsFrom(base):
            raise ValueError(f"Missing {base} {relative}:{key}")
        info = {"class": str(obj.ClassName())}
        if dimension is not None:
            if obj.GetDimension() != dimension or obj.InheritsFrom("TH2Poly"):
                raise ValueError(f"Unsupported histogram in {relative}:{key}")
            info["dimension"] = dimension
        self.files[relative]["objects"][key] = info
        return obj


def line(values):
    return " ".join(repr(float(value)) for value in values)


def axis_spec(axis):
    bins = axis.GetNbins()
    xbins = axis.GetXbins()
    if xbins.GetSize():
        edges = [float(xbins.At(i)) for i in range(xbins.GetSize())]
        uniform = 0
    else:
        edges = [float(axis.GetBinLowEdge(i)) for i in range(1, bins + 2)]
        edges[0], edges[-1] = float(axis.GetXmin()), float(axis.GetXmax())
        uniform = 1
    if len(edges) != bins + 1 or not all(math.isfinite(x) for x in edges):
        raise ValueError("Inconsistent ROOT axis")
    if not all(a < b for a, b in zip(edges, edges[1:])):
        raise ValueError("Unordered ROOT axis")
    return f"axis {uniform} {len(edges)} {line(edges)}"


def histogram1(obj):
    return [axis_spec(obj.GetXaxis()), line(obj.GetBinContent(i) for i in range(obj.GetNbinsX() + 2))]


def histogram2(obj):
    nx, ny = obj.GetNbinsX(), obj.GetNbinsY()
    count = (nx + 2) * (ny + 2)
    if obj.GetNcells() != count:
        raise ValueError("Unsupported ROOT TH2 storage")
    return [axis_spec(obj.GetXaxis()), axis_spec(obj.GetYaxis()),
            line(obj.GetBinContent(i) for i in range(count))]


def export_jet(sources):
    raw_config = f"JetCalibTools/{JET_AREA}CalibrationConfigs/{JET_CONFIG}"
    config = sources.root.TEnv(str(sources.resolve(raw_config)))
    settings = {
        key: str(config.GetValue(key, "")) for key in (
            "SmearType", "SmearingCorrectionHistType", "SmearingCorrectionInterpType",
            "SmearingCorrectionFile", "SmearingHistNameResolutionData", "SmearingHistNameResolutionMC",
        )
    }
    check_fields(settings, {
        "SmearType": "FourVec", "SmearingCorrectionHistType": "PtAbsEta",
        "SmearingCorrectionInterpType": "OnlyX",
        "SmearingCorrectionFile": "JetCalibTools/CalibrationFactors/JER_Nominal_Aug2018.root",
        "SmearingHistNameResolutionData": "JER_Nominal_data_AntiKt4EMPFlow",
        "SmearingHistNameResolutionMC": "JER_Nominal_MC_AntiKt4EMPFlow",
    })
    sources.files[raw_config]["settings"] = settings
    original = settings["SmearingCorrectionFile"]
    relative = original[:14] + JET_AREA + original[14:]
    rows = ["METSignificanceJetV1"]
    with sources.root_file(relative) as source:
        for tag, setting in (("data", "SmearingHistNameResolutionData"), ("mc", "SmearingHistNameResolutionMC")):
            key = settings[setting]
            h = sources.object(source, relative, key, "TH1", 2)
            x, y = h.GetXaxis(), h.GetYaxis()
            limits = [a.GetBinLowEdge(i) for a in (x, y) for i in (1, a.GetNbins() + 1)]
            if x.CanExtend() or y.CanExtend():
                raise ValueError(f"Dynamic axes are unsupported for {key}")
            rows += [tag + " " + line(limits), axis_spec(y)]
            info = sources.files[relative]["objects"][key]
            info["x_axis_range"] = [x.GetFirst(), x.GetLast(), bool(x.TestBit(sources.root.TAxis.kAxisRange))]
            info["projection_y_bins"] = y.GetNbins() + 1
            for bin_y in range(y.GetNbins() + 1):
                proj = h.ProjectionX(f"metsig_export_{tag}_{bin_y}", bin_y, bin_y)
                if not proj or proj.GetXaxis().CanExtend():
                    raise ValueError(f"Unsupported x projection for {key} y={bin_y}")
                proj.SetDirectory(0)
                sources.root.SetOwnership(proj, True)
                rows += ["projection", *histogram1(proj),
                         line(proj.GetBinCenter(i) for i in range(proj.GetNbinsX() + 2))]
    return rows


def export_egamma(sources):
    rows = ["METSignificanceEgammaV1"]
    suffix = "run2_release21_es2017_R21_v1.root"
    base = "ElectronPhotonFourMomentumCorrection/v20"
    count = None
    for kind, stem in (("electron", "electron"), ("unconverted", "recoUnconv"), ("converted", "recoConv")):
        relative = f"{base}/resolutionFit_{stem}_{suffix}"
        with sources.root_file(relative) as source:
            if kind == "electron":
                edge_hist = sources.object(source, relative, "hsamplingG", "TH1", 1)
                axis = edge_hist.GetXaxis()
                axis_spec(axis)  # Validate bounds before exporting this variable axis verbatim.
                bins = axis.GetXbins()
                if bins.GetSize() < 2:
                    raise ValueError("Stock eg_resolution requires electron hsamplingG variable edges")
                count = bins.GetSize() - 1
                rows.append(f"eta_edges {count + 1} " + line(bins.At(i) for i in range(count + 1)))
            rows.append(kind)
            for key in ("hsampling90", "hnoise90", "hconst90"):
                hist = sources.object(source, relative, key, "TH1", 1)
                if hist.GetNbinsX() != count:
                    raise ValueError(f"SigmaEff90 table has a different number of eta bins: {relative}:{key}")
                rows.append(line(hist.GetBinContent(i) for i in range(1, count + 1)))
    relative = "ElectronPhotonFourMomentumCorrection/v38/egammaEnergyCorrectionData.root"
    with sources.root_file(relative) as source:
        key = f"Resolution/{EGAMMA_MODEL}/ctZee_errStat"
        rows += ["photon_ct", *histogram1(sources.object(source, relative, key, "TH1", 1))]
    return rows


def export_muon(sources):
    prefix = f"MuonMomentumCorrections/{MUON_RELEASE}/ScaleAndSmear/Data22"
    rows = ["METSignificanceMuonData22V1"]
    relative = f"{prefix}/MC_values_muons.root"
    with sources.root_file(relative) as source:
        for tag, objects in (
            ("id_expected", ("r1_ID_MC", "r2_ID_MC", "r2_ID_TAN2_MC")),
            ("ms_expected", ("r0_MS_MC", "r1_MS_MC", "r2_MS_MC")),
        ):
            rows.append(tag)
            for key in objects:
                rows += histogram2(sources.object(source, relative, key, "TH2", 2))
    for tag, track, objects in (
        ("id_nominal", "ID", ("r1_nom", "r2_nom")),
        ("ms_nominal", "ME", ("r0_nom", "r1_nom", "r2_nom")),
    ):
        relative = f"{prefix}/{track}.root"
        rows.append(tag)
        with sources.root_file(relative) as source:
            for key in objects:
                rows += histogram2(sources.object(source, relative, key, "TH2", 2))
    return rows


def export_tau(sources):
    rows = ["METSignificanceTauV1"]
    with sources.root_file(TAU_FILE) as source:
        for mode in MODES:
            for eta_index in range(5):
                # Both objects are required by stock initialization. Only FinalCalib is evaluated.
                calo_key = f"tauRec/Graph_from_ResolutionEt_tauRec_{mode}_{eta_index}"
                sources.object(source, TAU_FILE, calo_key, "TGraph")
                key = f"FinalCalib/Graph_from_ResolutionEt_FinalCalib_{mode}_{eta_index}"
                graph = sources.object(source, TAU_FILE, key, "TGraph")
                n = graph.GetN()
                x = [float(graph.GetPointX(i)) for i in range(n)]
                y = [float(graph.GetPointY(i)) for i in range(n)]
                sorted_bit = bool(graph.TestBit(sources.root.TGraph.kIsSortedX))
                if not x or not all(math.isfinite(value) for value in x):
                    raise ValueError(f"Invalid TGraph abscissa in {key}")
                if sorted_bit and any(a > b for a, b in zip(x, x[1:])):
                    raise ValueError(f"Invalid TGraph sorted flag in {key}")
                maximum = max(x)
                rows += [f"graph {int(sorted_bit)} {maximum!r} {n}", line(x), line(y)]
                sources.files[TAU_FILE]["objects"][key].update(
                    sorted_x_flag=sorted_bit, points=n, max_pt_mev=maximum,
                )
    return rows


def extension(name, variables, table=None, outputs=None):
    data = {
        "nodetype": "extension",
        "type": f"atlas.met.{name}",
        "version": 1,
        "config": {"table": table} if table is not None else {},
    }
    output = {"name": "value", "type": "real"}
    if outputs:
        output.update(name="components", description=", ".join(outputs), shape=[len(outputs)])
        data["output_shape"] = [len(outputs)]
    return {
        "name": name, "version": 1, "output": output, "data": data,
        "inputs": [{"name": variable, "type": kind} for variable, kind in variables],
    }


def definitions(data):
    real, integer = "real", "int"
    return {"schema_version": 3, "corrections": [
        extension("jet_resolution", [("pt_mev", real), ("eta", real)], data["jet"],
                  ["data", "mc", "relative_pt_binary32"]),
        extension("egamma_resolution", [
            ("kind", integer), ("energy_mev", real), ("cluster_eta", real),
            ("eta_calo", real), ("has_eta_calo", integer), ("with_constant_term", integer),
        ], data["egamma"]),
        extension("muon_resolution", [
            ("detector", integer), ("run_number", integer), ("is_simulation", integer),
            ("eta", real), ("phi", real), ("cb_pt_gev", real), ("id_pt_gev", real), ("ms_pt_gev", real),
        ], data["muon"]),
        extension("tau_resolution", [
            ("decay_mode", integer), ("has_decay_mode", integer), ("tes_eta", real),
            ("final_eta", real), ("final_pt_mev", real),
        ], data["tau"]),
        extension("significance_finish", [(name, real) for name in (
            "varl_gev2", "vart_gev2", "cov_gev2", "met_gev", "sumet_gev", "soft_sumet_gev",
            "previous_squared", "previous_rho",
        )], outputs=["significance", "directional", "met_over_sqrt_sumet", "met_over_sqrt_ht",
                    "new_squared", "new_rho"]),
    ]}


def main(
    reference_config: Annotated[Path, typer.Option(exists=True, dir_okay=False)],
    output_dir: Annotated[Path, typer.Option()],
    algorithm: Annotated[str, typer.Option()] = "seq1_MetSignificanceAlg_AnaMET",
):
    """Extract the data and write a fresh CorrectionLib v4 plugin bundle."""
    import ROOT

    release = os.environ.get("AnalysisBase_VERSION") or os.environ.get("AtlasVersion")
    if release != "25.2.117" or ROOT.gROOT.GetVersion().replace("/", ".") != "6.40.02":
        raise ValueError("Source the AnalysisBase 25.2.117 environment before exporting")
    if ROOT.gSystem.Load("libPathResolverDict") < 0:
        raise RuntimeError("Could not load the release PathResolver dictionary")
    config = json.loads(reference_config.read_bytes())
    validate_config(config, algorithm)
    if output_dir.exists():
        raise FileExistsError(output_dir)
    sources = Sources(ROOT)
    data = {name: "\n".join(export(sources)) + "\n" for name, export in (
        ("jet", export_jet), ("egamma", export_egamma), ("muon", export_muon), ("tau", export_tau),
    )}
    output_dir.mkdir(parents=True, exist_ok=False)
    files = {}
    for name, encoded in data.items():
        path = output_dir / f"{name}.table"
        path.write_text(encoded, encoding="ascii")
        files[path.name] = sha256(path)
    payload = output_dir / "corrections.json"
    payload.write_text(json.dumps(definitions(data), separators=(",", ":")) + "\n", encoding="utf-8")
    files[payload.name] = sha256(payload)
    provenance = {
        "release": release, "root": ROOT.gROOT.GetVersion(), "athena": ATHENA,
        "correctionlib": CORE, "plugin_api": 4,
        "config": str(reference_config.resolve()), "config_sha256": sha256(reference_config),
        "algorithm": algorithm, "calibrations": sources.files, "outputs_sha256": files,
    }
    (output_dir / "provenance.json").write_text(
        json.dumps(provenance, indent=2, allow_nan=False) + "\n", encoding="utf-8",
    )


if __name__ == "__main__":
    typer.run(main)
