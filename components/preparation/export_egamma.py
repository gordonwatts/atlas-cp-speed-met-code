#!/usr/bin/env python3
"""Export release egamma leakage formulas, graphs and the pointing profile."""
import json
import math
import os
import struct
from pathlib import Path
from typing import Annotated

import typer

from egamma_schema import INTEGER_INPUTS, OUTPUTS, REAL_INPUTS
from export_tau import axis_line, check, digest, numbers

FILES = {
    "electron": ("IsolationCorrections/v6/isolation_ptcorrections_rel22_mc20.root",
                 "168cb3107ba39a0517b7182c343e484d3be65164ac609e9927f88efea59298dc"),
    "photon": ("IsolationCorrections/v5/isolation_ptcorrections_rel20_2.root",
               "5ce6d839cf194ecbc79dc79cb5f3a1172e719d0c65414319f541d9fcfbd2a0a2"),
    "profile": ("ElectronIsolationSelection/v1/CaloDeltaRZ.root",
                "16e94100e840af756df8a7fab0f0c21828bcb67e79aa033c2d9083b5415c04a2"),
}
ROOT_EXPRESSIONS = {
    "electron": "pow(10,[p0])*pow(x,[p1])*pow(x,[p2]*TMath::Log10(x))",
    "photon": "([p0]+[p1]*x+[p2]*pow(x,2))",
}
NATIVE_EXPRESSIONS = {
    "electron": "pow(10,[0])*pow(x,[1])*pow(x,[2]*log10(x))",
    "photon": "[0]+[1]*x+[2]*(x*x)",
}
ETA = ("0.0_0.1", "0.1_0.6", "0.6_0.8", "0.8_1.15", "1.15_1.37", "1.37_1.52",
       "1.52_1.81", "1.81_2.01", "2.01_2.37", "2.37_2.47")
PART = ("unconverted", "converted_ok", "converted_trouble", "author_1_electron", "author_16_electron")
CONE = (20, 30, 40)

def validate_config(path):
    config = json.loads(path.read_bytes())
    eta = "selectEta,as_bits"
    oq = "goodOQ,as_bits"
    clean = "isClean,as_bits"
    pt = "selectPt,as_bits"
    track = "trackSelection_loose,as_bits"
    final = {}
    for particle in ("Electron", "Photon"):
        suffix = f"_Ana{particle}s"
        photon = particle == "Photon"
        leading = "&&".join((eta, oq, clean) if photon else (eta, oq))
        pre = leading + "&&" + pt
        wp = "tight" if photon else "loose"
        wp_suffix = suffix + "_" + wp
        pid = "selectEM_tight,as_char" if photon else "selectLikelihood_loose,as_char"
        iso = f"isolated_{wp},as_char"

        def alg(name, expected, tool="selectionTool", tool_expected=None):
            props = config["seq1_" + name]["Properties"]
            check(props, expected)
            if tool_expected is not None:
                check(props["Tools"][tool]["Properties"], tool_expected)

        common = {
            "minEta": 0, "maxPt": 0, "minRapidity": 0, "maxRapidity": 0,
            "etaGapLow": 0, "etaGapHigh": 0, "useConstituentMomentum": 0, "useDressedProperties": 0,
        }
        pt_conf = common | {"minPt": 10000, "maxEta": 0, "useClusterEta": 0}
        alg(particle+"EtaCutAlg"+suffix, {"preselection": "", "selectionDecoration": eta},
            tool_expected=common | {"minPt": 0, "maxEta": 2.37 if photon else 2.47, "useClusterEta": 1})
        alg(particle+"ObjectQualityAlg"+suffix, {"preselection": eta, "selectionDecoration": oq},
            tool_expected={"Mask": 34214 if photon else 1446})
        tool = config["seq1_"+particle+"ObjectQualityAlg"+suffix]["Properties"]["Tools"]["selectionTool"]["Properties"]
        if "AsgDeadHVCellRemovalTool/" not in tool["deadHVTool"]:
            raise ValueError("Unsupported egamma dead-HV tool")
        if photon:
            alg("PhotonCleaningAlg"+suffix, {"preselection": eta+"&&"+oq, "selectionDecoration": clean},
                tool_expected={"selectionFlags": "['DFCommonPhotonsCleaning']", "invertFlags": "[0]"})
            alg("PhotonOriginCorrectionAlg"+suffix, {"preselection": leading})
        alg(particle+"CalibrationAndSmearingAlg"+suffix, {"preselection": leading, "skipNominal": True})
        alg(particle+"PtCutAlg"+suffix, {"preselection": leading, "selectionDecoration": pt},
            tool_expected=pt_conf)
        alg(particle+"IsolationCorrectionAlg"+suffix, {"preselection": pre},
            tool="isolationCorrectionTool", tool_expected={
                "ToolVer": "REL21" if photon else "REL22", "CorrFile": FILES[particle.lower()][0],
                "AFII_corr": 0, "Apply_SC_leakcorr": 0, "Apply_ddshifts": 0, "Apply_etaEDParPU_correction": 0,
                "Apply_etaEDPar_mc_correction": 0, "Correct_etcone": 0, "FixTimingIssueInCore": 1,
                "ForcePartType": 0, "IsMC": 1, "LogLogFitForLeakage": 0, "Trouble_categories": 1,
                "EventInfoKey": "EventInfo",
            })
        if not photon:
            alg("ElectronTrackSelectionAlg"+wp_suffix, {
                "preselection": pre, "selectionDecoration": track, "maxD0Significance": 5,
                "maxDeltaZ0SinTheta": 0.5,
            })
        selection_pre = pre if photon else pre+"&&"+track
        id_alg = "PhotonIsEMSelectorAlg" if photon else "ElectronLikelihoodAlg"
        id_input = "DFCommonPhotonsIsEMTight" if photon else "DFCommonElectronsLHLooseBL"
        alg(id_alg+wp_suffix, {"preselection": selection_pre, "selectionDecoration": pid},
            tool_expected={"selectionFlags": f"['{id_input}']", "invertFlags": "[0]"})
        alg(particle+"IsolationSelectionAlg"+wp_suffix,
            {"preselection": selection_pre+"&&"+pid, "selectionDecoration": iso},
            tool_expected={"IsoDecSuffix": "", "ElectronWP": "Undefined" if photon else "Tight_VarRad",
                           "PhotonWP": "FixedCutTight" if photon else "Undefined",
                           "ElectronWPVec": "[]", "PhotonWPVec": "[]", "MuonWPVec": "[]"})
        alg("PtEtaSelectionAlg"+suffix+"_", {"preselection": pre, "selectionDecoration": "selectPtEta"},
            tool_expected=pt_conf)
        final[particle.lower() + "sSelection"] = "&&".join((selection_pre, pid, iso, "selectPtEta"))
    check(config["seq1_MetMakerAlg_AnaMET"]["Properties"], final)

def node(name, config, *children, size=0):
    result = {"nodetype": "extension", "type": "atlas.met.preparation." + name,
              "version": 1, "config": config}
    if children:
        result["children"] = list(children)
    if size:
        result["output_shape"] = [size]
    return result

def correction(name, reals, ints, data, description="", size=0):
    output = {"name": "value", "type": "real", "description": description}
    if size:
        output["shape"] = [size]
    return {"name": name, "version": 1, "inputs": [
        *({"name": key, "type": "real"} for key in reals),
        *({"name": key, "type": "int"} for key in ints),
    ], "output": output, "data": data}

def u32(value):
    return struct.unpack("=I", struct.pack("=f", value))[0]

def f32(value):
    return struct.unpack("=f", struct.pack("=f", value))[0]

def probe(probes, kind, key, x, y, result):
    probes.append(f"{kind} {key} {x!r} {y!r} {u32(result)}")

def category(payload):
    return {"nodetype": "category", "input": "leak_key", "content": payload}

def extract_leakage(root, f, particle, probes):
    fit, low = [], []
    for ir, cone in enumerate(CONE):
        for j, eta in enumerate(ETA):
            for k in range(3) if particle == "photon" else range(3, 5):
                if k < 3 and j in (5, 9):
                    continue
                coord = f"{cone}_eta_{eta}_{PART[k]}"
                key = 50 * ir + 10 * k + j
                fit_key = "mean_f_topoetcone" + coord
                src = f.Get(fit_key)
                if not src or src.ClassName() != "TF1" or src.GetNpar() != 3 or (
                    str(src.GetExpFormula()).replace(" ", "") != ROOT_EXPRESSIONS[particle]
                ):
                    raise ValueError(f"Unsupported release leakage TF1 {fit_key}")
                parameters = [float(src.GetParameter(i)) for i in range(3)]
                if any(not math.isfinite(p) for p in parameters):
                    raise ValueError(f"Nonfinite TF1 parameter {fit_key}")
                fit.append({"key": key, "value": {
                    "nodetype": "formula", "expression": NATIVE_EXPRESSIONS[particle],
                    "parser": "TFormula", "variables": ["leak_et"], "parameters": parameters,
                }})
                for x in (500., 9999., 10000., 20000., 39251.37, 249999.999, 250000., 900000.):
                    probe(probes, particle+"_leakage_fit", key, x, 0., src.Eval(x))
                if k not in (2, 4):
                    continue
                graph_key = "mean_g_histo_topoetcone" + coord
                graph = f.Get(graph_key)
                if not graph or not graph.InheritsFrom("TGraph") or graph.GetN() < 1:
                    raise ValueError(f"Missing release leakage TGraph {graph_key}")
                n = graph.GetN()
                xs = [float(graph.GetPointX(i)) for i in range(n)]
                ys = [float(graph.GetPointY(i)) for i in range(n)]
                sorted_bit = int(graph.TestBit(root.TGraph.kIsSortedX))
                if any(not math.isfinite(x) for x in xs+ys):
                    raise ValueError(f"Nonfinite source graph point {graph_key}")
                encoded = f"MetLeakageGraphV1 {sorted_bit} {n} {numbers(xs)} {numbers(ys)}"
                low.append({"key": key, "value": node("leakage_graph", {"points": encoded})})
                for x in [min(xs)-1000., max(xs)+1000., *xs, *[(a+b)/2 for a,b in zip(xs, xs[1:])]]:
                    probe(probes, particle+"_leakage_low", key, x, 0., graph.Eval(x))
    return category(fit), category(low)

def extract_pointing(f, probes):
    hist = f.Get("hData")
    if not hist or hist.ClassName() != "TProfile2D" or hist.GetDimension() != 2:
        raise ValueError("Expected release hData TProfile2D")
    axes = hist.GetXaxis(), hist.GetYaxis()
    metadata = ["MetCaloDepthV1"]
    for axis in axes:
        n = axis.GetNbins()
        metadata.extend((axis_line(axis),
                         numbers(axis.GetBinCenter(i) for i in range(n+2)),
                         numbers(axis.GetBinUpEdge(i) for i in range(1, n+1)),
                         numbers(axis.GetBinWidth(i) for i in range(1, n+1))))
    nx, ny = (axis.GetNbins() for axis in axes)
    regular_content = [float(hist.GetBinContent(hist.GetBin(i, j))) for i in range(1, nx+1)
                       for j in range(1, ny+1)]
    table = {
        "nodetype": "multibinning", "inputs": ["rz_x", "rz_y"], "flow": "error",
        "edges": [[float(axis.GetBinLowEdge(i)) for i in range(1, axis.GetNbins()+2)] for axis in axes],
        "content": regular_content,
    }
    def stock_limit(axis, value):
        requested = axis.FindFixBin(value)
        if requested < 1:
            return f32(axis.GetBinLowEdge(1) + f32(1.e-6))
        if requested > axis.GetNbins():
            return f32(axis.GetBinUpEdge(axis.GetNbins()) - f32(1.e-6))
        return value
    for dimension, axis in enumerate(axes):
        size = axis.GetNbins()
        centres = [axis.GetBinCenter(k) for k in range(size+2)]
        edges = [axis.GetBinLowEdge(k) for k in range(1, size+2)]
        for test in centres + edges:
            coords = [f32(0.17), f32(-0.97)]
            coords[dimension] = f32(test)
            limited = [stock_limit(s, v) for s, v in zip(axes, coords)]
            probe(probes, "egamma_pointing_depth", 0, *coords, hist.Interpolate(*limited))
    return node("pointing_profile", {"depth_axes": "\n".join(metadata)}, table)

def export(
    reference_config: Annotated[Path, typer.Option(exists=True, dir_okay=False, readable=True)],
    output_dir: Annotated[Path, typer.Option(file_okay=False)],
):
    """Extract offline egamma calibration objects and independent ROOT probes."""
    validate_config(reference_config)
    import ROOT
    ROOT.PyConfig.IgnoreCommandLineOptions = True
    if ROOT.gSystem.Load("libPathResolver") < 0:
        raise RuntimeError("Load the AnalysisBase release before exporting")
    resolved = {}
    for kind, (asset, expected_hash) in FILES.items():
        value = str(ROOT.PathResolverFindCalibFile(asset))
        if not value:
            raise FileNotFoundError(asset)
        path = Path(value).resolve(strict=True)
        if digest(path) != expected_hash:
            raise ValueError(f"Changed source calibration hash for {asset}")
        resolved[kind] = path

    probes, corrections = [], []
    source = ROOT.TFile.Open(str(resolved["profile"]), "READ")
    if not source or source.IsZombie():
        raise RuntimeError("Cannot open shower-depth profile")
    try:
        profile = extract_pointing(source, probes)
    finally:
        source.Close()
    corrections.append(correction("egamma_pointing_depth", ("rz_x", "rz_y"), (), profile))
    for particle in ("electron", "photon"):
        source = ROOT.TFile.Open(str(resolved[particle]), "READ")
        if not source or source.IsZombie():
            raise RuntimeError(f"Cannot open {particle} leakage source")
        try:
            fit, low = extract_leakage(ROOT, source, particle, probes)
        finally:
            source.Close()
        corrections.extend((
            correction(particle+"_leakage_fit", ("leak_et",), ("leak_key",), fit),
            correction(particle+"_leakage_low", ("leak_et",), ("leak_key",), low),
            correction(particle+"_prepare", REAL_INPUTS, INTEGER_INPUTS,
                       node("egamma", {"photon": int(particle == "photon"), "recipe": "25.2.117-fixed-fullsim"},
                            fit, low, profile, size=len(OUTPUTS)), ", ".join(OUTPUTS), len(OUTPUTS)),
        ))

    output_dir.mkdir(parents=True, exist_ok=True)
    payload = output_dir / "egamma.json"
    temporary = payload.with_suffix(".tmp")
    temporary.write_text(json.dumps({
        "schema_version": 2, "description": "AnalysisBase 25.2.117 fixed nominal MET egamma preparation",
        "corrections": corrections,
    }, allow_nan=False, separators=(",", ":")) + "\n")
    os.replace(temporary, payload)
    table = output_dir / "egamma-probes.tsv"
    table.write_text("correction key x y result_float_bits\n" + "\n".join(probes) + "\n")
    manifest = {
        "athena": "6ce3f7ae7c41dd9eaec170c3c349df15faef7e5c", "root_version": ROOT.gROOT.GetVersion(),
        "calibrations": {kind: {"path": str(path), "sha256": digest(path)} for kind, path in resolved.items()},
        "configuration": {"path": str(reference_config), "sha256": digest(reference_config)},
        "payload_sha256": digest(payload), "probe_sha256": digest(table), "probes": len(probes),
    }
    (output_dir / "egamma-manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    typer.echo(f"Wrote {payload} and {len(probes)} probes")

if __name__ == "__main__":
    typer.run(export)
