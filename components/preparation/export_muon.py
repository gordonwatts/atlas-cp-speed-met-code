#!/usr/bin/env python3
"""Export the initialized Run 3 nominal MET muon selections."""
import hashlib
import json
from pathlib import Path
from typing import Annotated

import typer

from export_tau import check

PREFIX = "seq1_"
ETA = "selectEta,as_bits"
PT = "selectPt,as_bits"
TRACK = "trackSelection_medium,as_bits"
MEDIUM = "good_muon_medium,as_char"
ISOLATED = "isolated_muon_medium,as_char"

def joined(*cuts):
    return "&&".join(cuts)

def export(
    reference_config: Annotated[Path, typer.Option(exists=True, dir_okay=False, readable=True)],
    output_dir: Annotated[Path, typer.Option(file_okay=False)],
):
    """Check the release recipe and write a muon payload with source hashes."""
    config = json.loads(reference_config.read_bytes())

    def alg(name, properties, tools=None):
        stage = config[PREFIX + name]["Properties"]
        check(stage, properties)
        for tool, expected in (tools or {}).items():
            check(stage["Tools"][tool]["Properties"], expected)

    common = {
        "useClusterEta": 0, "useConstituentMomentum": 0, "useDressedProperties": 0,
        "etaGapHigh": 0, "etaGapLow": 0, "minEta": 0, "maxPt": 0,
        "minRapidity": 0, "maxRapidity": 0,
    }
    alg("MuonEtaCutAlg_AnaMuons", {"preselection": "", "selectionDecoration": ETA},
        {"selectionTool": common | {"maxEta": 2.7, "minPt": 0}})
    alg("MuonCalibrationAndSmearingAlg_AnaMuons", {"skipNominal": True, "preselection": ETA})
    alg("MuonPtCutAlg_AnaMuons", {"preselection": ETA, "selectionDecoration": PT},
        {"selectionTool": common | {"maxEta": 0, "minPt": 3000}})
    alg("MuonTrackSelectionAlg_AnaMuons_medium", {
        "preselection": joined(ETA, PT), "selectionDecoration": TRACK,
        "maxD0Significance": 3, "maxDeltaZ0SinTheta": 0.5,
    })
    alg("MuonSelectionAlg_AnaMuons_medium", {
        "preselection": joined(ETA, PT, TRACK), "selectionDecoration": MEDIUM,
        "badMuonVetoDecoration": "is_bad_medium,as_char",
    }, {"selectionTool": {
        "MuQuality": 1, "MaxEta": 2.7, "IsRun3Geo": 1, "AllowSettingGeometryOnTheFly": 0,
        "RecalcPrecisionLayerswNSW": 1, "ExcludeNSWFromPrecisionLayers": 0, "TurnOffMomCorr": 0,
        "CaloScoreWP": 4, "UseCaloScore": 1, "UseAllAuthors": 1, "ToroidOff": 0,
        "DoBadMuonVetoMimic": 0, "DisablePtCuts": 0, "UseLRT": 0,
        "PixCutOff": 0, "SctCutOff": 0, "SiHolesCutOff": 0, "TrtCutOff": 1,
        "UseMVALowPt": 0, "UseSegmentTaggedLowPt": 0, "ExpertDevelopMode": 0,
        "CalculateTightNNScore": 0, "Use2stationMuonsHighPt": 1, "UseBEEBISInHighPtRun3": 1,
    }})
    alg("MuonIsolationAlg_AnaMuons_medium", {
        "preselection": joined(ETA, PT, TRACK, MEDIUM), "isolationDecoration": ISOLATED,
    }, {"isolationTool": {
        "MuonWP": "Loose_VarRad", "MuonWPVec": "[]", "IsoDecSuffix": "",
        "ElectronWP": "Undefined", "PhotonWP": "Undefined",
    }})
    alg("MetMakerAlg_AnaMET", {"muonsSelection": joined(ETA, PT, TRACK, MEDIUM, ISOLATED)})

    local = Path(__file__).resolve().parent
    graph = local / "muon_schema.json"
    raw = graph.read_bytes()
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "muon.json").write_bytes(raw)
    sources = [reference_config, graph, local / "muon.cxx", Path(__file__).resolve()]
    manifest = {
        "recipe": "AnalysisBase-25.2.117-Run3-Medium-Loose_VarRad-NOSYS",
        "athena": "6ce3f7ae7c41dd9eaec170c3c349df15faef7e5c",
        "payload_sha256": hashlib.sha256(raw).hexdigest(),
        "sources": {path.name: hashlib.sha256(path.read_bytes()).hexdigest() for path in sources},
        "isolation": {
            "track": "ptvarcone30_Nonprompt_All_MaxWeightTTVA_pt1000",
            "calo": "topoetcone20",
        },
    }
    (output_dir / "muon-manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    typer.echo(f"Wrote {output_dir / 'muon.json'}")

if __name__ == "__main__":
    typer.run(export)
