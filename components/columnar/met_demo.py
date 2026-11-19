# /// script
# requires-python = ">=3.14,<3.15"
# dependencies = ["atlas-met-columnar==0.1.0", "awkward==2.15.0", "hist", "numpy", "typer", "uproot==5.7.7"]
# [tool.uv.sources]
# atlas-met-columnar = {path = "wheels/atlas_met_columnar-0.1.0-cp314-cp314-linux_x86_64.whl"}
# ///
"""Run all five frozen MET states directly on original PHYSLITE."""
from __future__ import annotations

import hashlib
from importlib.metadata import version
import json
from pathlib import Path
import sys

try:
    import atlas_met_portable as installed_runtime
except ModuleNotFoundError as error:
    if error.name != "atlas_met_portable":
        raise
    installed_runtime = None

import hist
import numpy as np
import typer
import uproot

from atlas_met_columnar import ASSOC, ASSOC_FIELDS, RAW_BRANCHES, MemberwiseJets, Native, pack_raw, read_associations
from atlas_met_prepare import ADDITIONAL, prepare_chunk

UUID = "79ee6652-311f-11f1-aa37-0445010abeef"
STATES = ["NOSYS", "MET_SoftTrk_ResoPara", "MET_SoftTrk_ResoPerp",
          "MET_SoftTrk_Scale__1down", "MET_SoftTrk_Scale__1up"]
TERM_NAMES = ["RefEle", "RefGamma", "RefTau", "Muons", "MuonEloss", "RefJet", "PVSoftTrk", "Final"]


def book(labels):
    return hist.Hist(hist.axis.StrCategory(labels, overflow=False, name="state"),
                     hist.axis.Regular(100, 0, 500, name="met_GeV", underflow=True, overflow=True),
                     storage=hist.storage.Weight())


def signature(histogram):
    return hashlib.sha256(histogram.values(flow=True).tobytes()).hexdigest()


def save(histogram, labels, entries, output):
    sumw = histogram.values(flow=True)
    sumw2 = histogram.variances(flow=True)
    counts = sumw.astype("i8")
    np.testing.assert_array_equal(counts, sumw)
    np.testing.assert_array_equal(counts, sumw2)
    np.testing.assert_array_equal(counts.sum(axis=1), np.full(len(labels), entries))
    with output.open("xb") as file:
        np.savez(file, state_labels=np.asarray(labels), bin_edges_gev=histogram.axes[1].edges,
                 counts=counts, sumw=sumw, sumw2=sumw2)


def read_owner_keys(runtime):
    metadata = json.loads((runtime/"payloads"/"EventFormat.json").read_text())
    if metadata["input_uuid"] != UUID:
        raise ValueError("Wrong original EventFormat owner inventory")
    owners = {name: item["key"] for name,item in metadata["top_level_eventformat_elements"].items()}
    expected = {"AnalysisJets":452794649,"AnalysisElectrons":956497600,"AnalysisTauJets":518718875,
                "AnalysisMuons":980095599,"AnalysisPhotons":902907695}
    if any(owners.get(name)!=value for name,value in expected.items()):
        raise ValueError("EventFormat owner keys disagree with the memberwise association contract")
    return owners


def save_runtime_record(output, runtime, input_layout):
    paths = set()
    for row in Path("/proc/self/maps").read_text().splitlines():
        fields = row.split(maxsplit=5)
        if len(fields) == 6 and fields[5].startswith("/"):
            paths.add(fields[5])
    shared_objects = sorted(path for path in paths if ".so" in Path(path).name)
    forbidden = tuple(path for path in shared_objects if "/cvmfs/" in path or
                      Path(path).name.startswith(("libCore.so","libRIO.so","libTree.so","libHist.so",
                                                  "libMathCore.so","libxAOD","libAthena")))
    required = {"libatlas_met_columnar.so", "libatlas_met_maker_correction.so",
                "libatlas_met_maker_provider.so", "libmet_preparation_provider.so", "libcorrectionlib.so"}
    basenames = {Path(path).name for path in shared_objects}
    if forbidden or not required <= basenames:
        raise RuntimeError(f"Unexpected live libraries: forbidden={forbidden}, missing={required-basenames}")
    payloads = {name.name: hashlib.sha256(name.read_bytes()).hexdigest()
                for name in sorted((runtime/"payloads").glob("*.json"))}
    dependencies = {name:version(name) for name in ("awkward","uproot","numpy","hist","typer")}
    with output.open("x") as handle:
        json.dump({"python":sys.version,"dependencies":dependencies,"input_uuid":UUID,
                   "input_layout":input_layout,"payload_sha256":payloads,
                   "mapped_shared_objects":shared_objects},handle,indent=2)
        handle.write("\n")


def main(
    input_file: Path, output: Path,
    runtime: Path | None = None,
    repeat: int = typer.Option(1, "-n", "--repeat", min=1),
    stop: int | None = typer.Option(None, min=0),
    chunk_size: int = typer.Option(2000, min=1),
    nominal_only: bool = False,
    correlated_diagnostic: bool = False,
    trace: Path | None = None,
    observables: Path | None = None,
    runtime_record: Path | None = None,
):
    if sys.version_info[:2] != (3, 14):
        raise ValueError("This demonstration requires CPython 3.14")
    for path in (output,observables,trace,runtime_record):
        if path is not None and path.exists():
            raise FileExistsError(path)
    if repeat != 1 and (trace or observables):
        raise ValueError("Per-event diagnostics require one input pass")
    if correlated_diagnostic and nominal_only:
        raise ValueError("Correlated resolution is a separate diagnostic mode")
    if runtime is None:
        if installed_runtime is None:
            raise ValueError("Install the paired wheel or supply the source-checkpoint --runtime directory")
        runtime = Path(installed_runtime.__file__).resolve().parent
    runtime = runtime.resolve(strict=True)
    native_path = runtime/"native"
    payloads = runtime/"payloads"
    keys = read_owner_keys(runtime)
    bridge = Native(native_path/"libatlas_met_columnar.so")
    # This C++ ABI loads the installed pinned CorrectionLib core, the shared API-v4
    # registry and named corrections. Each evaluate below enters CorrectionLib.
    context = bridge.open(native_path/"libatlas_met_maker_provider.so",
                          payloads/"maker.json",payloads/"nnjvt.json")
    labels = ["NOSYS"] if nominal_only else STATES.copy()
    if correlated_diagnostic:
        labels = ["NOSYS","MET_SoftTrk_ResoCorr"]
    if trace:
        trace = trace.resolve()
        with trace.open("x"):
            pass
    obs_terms, obs_sources, obs_sig = [], [], []
    obs_ids = {"event.number": [], "event.run": [], "event.dsid": []}
    expected_digest = None
    input_layout = None
    try:
        bridge.configure(context,payloads/"soft.json",payloads/"significance.json",
                         native_path/"libmet_preparation_provider.so",payloads/"preparation.json")
        for _ in range(repeat):
            bridge.lib.met_columnar_reset(context)
            histogram = book(labels)
            entries = 0
            with uproot.open(input_file,handler=uproot.source.file.MemmapSource,
                             num_workers=1,num_fallback_workers=1) as original:
                if str(original.file.uuid) != UUID:
                    raise ValueError("This payload and the event-local link inventory are pinned to the original file")
                tree = original["CollectionTree"]
                end = tree.num_entries if stop is None else min(stop,tree.num_entries)
                fields = list(dict.fromkeys([*RAW_BRANCHES,*(f for f in ADDITIONAL if f in tree)]))
                link_reader = MemberwiseJets(bridge,tree[ASSOC][ASSOC+"jetLink"])
                if runtime_record is not None:
                    input_layout = {"tree":"CollectionTree","loaded_branches":fields,
                                    "supported_but_absent":[f for f in ADDITIONAL if f not in tree],
                                    "nested_associations":[ASSOC+f for f in ASSOC_FIELDS],
                                    "memberwise_basket":[ASSOC+"jetLink"],"split_leaves":{}}
                for first in range(0,end,chunk_size):
                    last = min(first+chunk_size,end)
                    raw = tree.arrays(fields,entry_start=first,entry_stop=last,library="ak",how=dict)
                    if runtime_record is not None and first == 0:
                        input_layout["split_leaves"] = {name: list(children) for name,children in raw.items()
                                                         if isinstance(children,dict)}
                    assoc = read_associations(tree,first,last)
                    buffers = pack_raw(raw,assoc,link_reader.read(first,last))
                    # Nominal corrections are explicitly selected by their declared names
                    # in the raw preparation projections and evaluated on each typed batch.
                    evaluate = lambda name, inputs: bridge.prepare(context,name,inputs)
                    prepare_chunk(buffers,raw,keys,evaluate,capture_preparation=trace is not None)
                    # The native chunk calls shared maker, correct_pvsofttrk, build_final,
                    # and met_significance_event in event order, with the explicit labels.
                    terms,sources,significance = bridge.full(context,buffers,first,last-first,labels,trace)
                    for j,label in enumerate(labels):
                        gev = terms[:,j,7,3].astype("f8")/1000.
                        if not np.isfinite(gev).all():
                            raise ValueError(f"Nonfinite Final MET in {label} at entries {first}:{last}")
                        histogram.fill(label,gev,threads=1)
                    if observables is not None:
                        obs_terms.append(terms)
                        obs_sources.append(sources)
                        obs_sig.append(significance)
                        for key, values in obs_ids.items():
                            values.append(buffers[key])
                    entries += last-first
            digest = signature(histogram)
            if expected_digest is not None and expected_digest != digest:
                raise RuntimeError("Repeated original file passes did not produce identical yields")
            expected_digest = digest
        if runtime_record is not None:
            save_runtime_record(runtime_record, runtime, input_layout)
    finally:
        bridge.lib.met_columnar_close(context)
    save(histogram,labels,entries,output)
    if observables is not None:
        with observables.open("xb") as handle:
            np.savez(handle, states=np.asarray(labels), term_names=np.asarray(TERM_NAMES),
                     term_fields=np.asarray(["mpx","mpy","sumet","met","phi"]),
                     terms=np.concatenate(obs_terms) if obs_terms else np.empty((0,len(labels),8,5),"f4"),
                     sources=np.concatenate(obs_sources) if obs_sources else np.empty((0,len(labels),8),"u8"),
                     significance=np.concatenate(obs_sig) if obs_sig else np.empty((0,len(labels),4),"f4"),
                     entry=np.arange(entries,dtype="u8"),
                     event_number=np.concatenate(obs_ids["event.number"]) if entries else np.empty(0,"u8"),
                     run_number=np.concatenate(obs_ids["event.run"]) if entries else np.empty(0,"u4"),
                     mc_channel_number=np.concatenate(obs_ids["event.dsid"]) if entries else np.empty(0,"u4"))
    print(json.dumps({"output":str(output),"input_uuid":UUID,"entries":entries,
                      "event_visits":repeat*entries,"repeats":repeat,"states":labels,
                      "histogram_sha256":expected_digest}))


if __name__ == "__main__":
    typer.run(main)
