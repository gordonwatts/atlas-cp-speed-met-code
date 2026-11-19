#!/usr/bin/env python3
"""Export and verify the frozen AnalysisBase 25.2.117 TST content."""

import hashlib
from itertools import product
import json
import math
import os
from pathlib import Path
import struct
from typing import Annotated

import typer

RELEASE = "25.2.117"
ATHENA = "6ce3f7ae7c41dd9eaec170c3c349df15faef7e5c"
CORRECTIONLIB = "3dcfa20ec01fb40c4d5fec58383785519dec1c00"
RECIPE = "AnalysisBase-25.2.117-PFlow-Tight-PVSoftTrk"
PREFIX = "METUtilities/R22_PreRecs/"
CONFIG = "TrackSoftTerms-pflow_Dec24.config"
HIST_FILE = PREFIX + "R22_PreRecs_Run2Run3_TSTSystematics_Dec24.root"
FIXED_CALIBRATION_HASH = {
    PREFIX + CONFIG: "368cb13d92c9d45fe7da4c807466a6a720303f22400a2f1b89e387a061d3ddcd",
    HIST_FILE: "d402bcdb1d68790d26176a924554567e28e9cace3cfb35d22d75dd3904bea987",
}
TABLES = ("shiftpara", "resopara", "resoperp")
AXES = ("pthard", "njet", "mu")
REAL_INPUTS = [(name, "real") for name in AXES]

app = typer.Typer(no_args_is_help=True)


def file_hash(path):
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def write_json(path, value):
    with path.open("x", encoding="utf-8") as handle:
        json.dump(value, handle, separators=(",", ":"), allow_nan=False)
        handle.write("\n")


def expect(props, expected):
    for name, value in expected.items():
        actual = str(props.get(name)).strip("'")
        if actual != value:
            raise ValueError(f"Unsupported setting {name}={actual!r}, expected {value!r}")


def validate_job(config, maker_name, builder_name):
    maker = config[maker_name]
    builder = config[builder_name]
    if maker["Type"] != "CP::MetMakerAlg" or builder["Type"] != "CP::MetBuilderAlg":
        raise ValueError("Expected stock MET maker and builder")
    expect(maker["Properties"], {"softTermKey": "PVSoftTrk"})
    expect(builder["Properties"], {"softTerm": "PVSoftTrk"})
    if "finalKey" in builder["Properties"]:
        expect(builder["Properties"], {"finalKey": "Final"})
    tool = maker["Properties"]["Tools"]["systematicsTool"]
    if tool["Type"] != "met::METSystematicsTool":
        raise ValueError("Expected stock METSystematicsTool")
    expect(tool["Properties"], {
        "ConfigPrefix": PREFIX, "ConfigSoftTrkFile": CONFIG,
        "ConfigSoftCaloFile": "", "ConfigJetTrkFile": "", "UseDevArea": "0",
        "TruthContainer": "MET_Truth", "EventInfo": "EventInfo",
    })


def calibration(root, logical):
    resolved = root.PathResolver.FindCalibFile(logical)
    if not resolved:
        raise RuntimeError(f"PathResolver could not resolve {logical}")
    path = Path(str(resolved)).resolve(strict=True)
    digest = file_hash(path)
    if digest != FIXED_CALIBRATION_HASH[logical]:
        raise ValueError(f"Unsupported source contents at {path}, SHA256={digest}")
    return path, {"logical_name": logical, "resolved": str(path), "sha256": digest,
                  "size": path.stat().st_size}


def axis_signature(axis):
    edges = [float(axis.GetBinLowEdge(i)) for i in range(1, axis.GetNbins() + 2)]
    variable = axis.GetXbins().GetSize() != 0
    if not all(math.isfinite(edge) for edge in edges) or any(a >= b for a, b in zip(edges, edges[1:])):
        raise RuntimeError("Malformed source axis")
    if axis.CanExtend():
        raise RuntimeError("Soft term table cannot preserve an extendable ROOT axis")
    return variable, tuple(edges)


def basis_table(axes, histogram):
    # The recorded source shares shift-table indices with the resolution tables.
    nx, ny, nz = (axis.GetNbins() for axis in axes)
    content = [
        float(histogram.GetBinContent(min(x, nx), y, z))
        for x, y, z in product(range(nx + 2), range(ny + 2), range(nz + 2))
    ]
    edges = [axis_signature(axis)[1] for axis in axes]
    return {
        "nodetype": "multibinning",
        "inputs": list(AXES),
        "edges": [["-inf", *edge, "inf"] for edge in edges],
        "content": content,
        "flow": "clamp",
    }


def correction(name, inputs, node, shape=None):
    output = {"name": "result", "type": "real"}
    if shape:
        output["shape"] = [shape]
    return {
        "name": name, "description": RECIPE, "version": 1,
        "inputs": [{"name": variable, "type": kind} for variable, kind in inputs],
        "output": output, "data": node,
    }


def extension(name, suffix, inputs, children=(), shape=None):
    node = {"nodetype": "extension", "type": f"atlas.met.{suffix}", "version": 1,
            "config": {"recipe": RECIPE}, "children": list(children)}
    if shape:
        node["output_shape"] = [shape]
    return correction(name, inputs, node, shape)


def stock_probes(axes, name, hist):
    counts = [axis.GetNbins() for axis in axes]
    default = [float(axis.GetBinCenter(1)) for axis in axes]

    def expected(coords, mode):
        regular = [axis.FindFixBin(val) for axis, val in zip(axes, coords)]
        regular[0] = min(regular[0], counts[0])
        return {"name": name, "coords": [float(value).hex() for value in coords],
                "bins": regular, "expected": float(hist.GetBinContent(*regular)).hex(), "mode": mode}

    for i, axis in enumerate(axes):
        points = [float(axis.GetBinCenter(b)) for b in range(counts[i] + 2)]
        for edge in axis_signature(axis)[1]:
            points.extend((math.nextafter(edge, -math.inf), edge, math.nextafter(edge, math.inf)))
        points.extend((-math.inf, math.inf, math.nan))
        for val in points:
            point = default.copy()
            point[i] = val
            yield expected(point, AXES[i])
    for bins in product(*(range(count + 2) for count in counts)):
        coords = [float(axis.GetBinCenter(index)) for axis, index in zip(axes, bins)]
        yield expected(coords, "all-bin-centers")


@app.command()
def export(
    reference_config: Annotated[Path, typer.Option(exists=True, dir_okay=False)],
    output_dir: Annotated[Path, typer.Option()],
    maker: Annotated[str, typer.Option()] = "seq1_MetMakerAlg_AnaMET",
    builder: Annotated[str, typer.Option()] = "seq1_MetBuilderAlg_AnaMET",
):
    """Resolve the frozen TEnv and TH3D files with stock PyROOT and PathResolver."""
    import ROOT

    if os.environ.get("AtlasVersion") != RELEASE or str(ROOT.gROOT.GetVersion()).replace("/", ".") != "6.40.02":
        raise RuntimeError("The AnalysisBase 25.2.117 release must be active")
    if ROOT.gSystem.Load("libPathResolverDict") < 0:
        raise RuntimeError("Could not load PathResolver")
    validate_job(json.loads(reference_config.read_bytes()), maker, builder)
    if output_dir.exists():
        raise FileExistsError(output_dir)
    config_path, config_source = calibration(ROOT, PREFIX + CONFIG)
    config = ROOT.TEnv(str(config_path))
    settings = {key: str(config.GetValue(key, "")) for key in
                ("Conf.InputFile", "Units", "PtHardSyst.sourcedir", "PtHardSyst.suffix")}
    expect(settings, {
        "Conf.InputFile": HIST_FILE, "Units": "", "PtHardSyst.sourcedir": "/", "PtHardSyst.suffix": "TST",
    })

    root_path, root_source = calibration(ROOT, HIST_FILE)
    source = ROOT.TFile.Open(str(root_path), "READ")
    if not source or source.IsZombie():
        raise RuntimeError(f"Cannot open original TH3D file {root_path}")
    try:
        objects = [f"{settings['PtHardSyst.sourcedir']}/{table}_{settings['PtHardSyst.suffix']}" for table in TABLES]
        histograms = [source.Get(key) for key in objects]
        if any(not hist or str(hist.ClassName()) != "TH3D" for hist in histograms):
            raise RuntimeError(f"Missing original TH3Ds {objects}")
        axes = [getattr(histograms[0], f"Get{axis}axis")() for axis in "XYZ"]
        signatures = [axis_signature(axis) for axis in axes]
        for hist in histograms[1:]:
            if signatures != [axis_signature(getattr(hist, f"Get{axis}axis")()) for axis in "XYZ"]:
                raise RuntimeError("Source resolution tables do not share shift table axis definitions")
        if [axis.GetNbins() for axis in axes] != [16, 1, 1]:
            raise RuntimeError("Unexpected source table dimensions")
        children = [basis_table(axes, hist) for hist in histograms]
        probes = [row for name, hist in zip(TABLES, histograms) for row in stock_probes(axes, name, hist)]
    finally:
        source.Close()

    standalone = {"schema_version": 2, "corrections": [
        correction(name, REAL_INPUTS, node) for name, node in zip(TABLES, children)
    ]}
    nodes = {"schema_version": 3, "corrections": [
        extension("correct_pvsofttrk", "soft_track", [("ticket", "int"), ("state", "string"), *REAL_INPUTS],
                  children=children),
        extension("build_final", "final_sum", [("ticket", "int")]),
        extension("met_phi", "met_phi", [("mpx", "real"), ("mpy", "real")], shape=2),
    ]}
    output_dir.mkdir(parents=True, exist_ok=False)
    files = {"corrections.json": nodes, "tables.json": standalone, "table-probes.json": probes}
    for name, obj in files.items():
        write_json(output_dir / name, obj)
    with (output_dir / "table-probes.tsv").open("x", encoding="ascii") as stream:
        for row in probes:
            stream.write("\t".join((row["name"], *row["coords"], row["expected"])) + "\n")
    write_json(output_dir / "provenance.json", {
        "analysisbase": RELEASE, "athena": ATHENA, "root": str(ROOT.gROOT.GetVersion()),
        "correctionlib_pin": CORRECTIONLIB, "plugin_api": 4, "recipe": RECIPE,
        "stock_config": str(reference_config.resolve()), "stock_config_sha256": file_hash(reference_config),
        "tenv": config_source, "tenv_settings": settings,
        "root_input": root_source, "histograms": objects, "axes": signatures,
        "units": "MeV", "probe_rows": len(probes),
        "outputs_sha256": {name: file_hash(output_dir / name) for name in (*files, "table-probes.tsv")},
    })


@app.command("check-tables")
def check_tables(
    bundle: Annotated[Path, typer.Option(exists=True, file_okay=False)],
    output: Annotated[Path, typer.Option()],
):
    """Compare stock table values exactly through the pinned native binding."""
    from correctionlib import _core

    payload = bundle / "tables.json"
    probes = json.loads((bundle / "table-probes.json").read_bytes())
    evaluators = {"file": _core.CorrectionSet.from_file(str(payload)),
                  "string": _core.CorrectionSet.from_string(payload.read_text())}
    issues = {}
    for loader, tables in evaluators.items():
        issues[loader] = []
        for check_id, probe in enumerate(probes):
            args = [float.fromhex(value) for value in probe["coords"]]
            expected = float.fromhex(probe["expected"])
            try:
                got = tables[probe["name"]].evaluate(*args)
            except (RuntimeError, ValueError) as exc:
                issues[loader].append({"id": check_id, "error": str(exc)})
                continue
            if struct.pack(">d", expected) != struct.pack(">d", got):
                issues[loader].append({"id": check_id, "actual": got.hex(), "probe": probe})
    core = Path(_core.__file__).resolve()
    write_json(output, {
        "expected_core_pin": CORRECTIONLIB, "core": str(core), "core_sha256": file_hash(core),
        "table_rows_per_loader": len(probes), "tables_sha256": file_hash(payload), "failures": issues,
    })
    if any(issues.values()):
        raise typer.Exit(1)


if __name__ == "__main__":
    app()
