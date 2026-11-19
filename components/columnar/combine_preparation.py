"""Merge the peer-exported native preparation and field-access definitions."""
import json
from pathlib import Path
import typer


def extension(key, kind, inputs, components=None, version=1):
    output = {"name":"value","type":"real"}
    data = {"nodetype":"extension","type":kind,"version":version,"config":{}}
    if components:
        output["shape"] = [components]
        data["output_shape"] = [components]
    return {"name":key,"version":version,"output":output,"data":data,
            "inputs":[{"name":name,"type":dtype} for name,dtype in inputs]}


def main(output: Path, inputs: list[Path]):
    corrections = [
        extension("met_prepared_p4","atlas.met.prepared_p4",
                  [(k,"real") for k in ("pt","eta","phi","mass")]+[("particle_kind","int")],3,version=2),
        extension("met_track_pt","atlas.met.track_pt",
                  [(k,"real") for k in ("qop","theta","phi")]+[("valid","int")]),
    ]
    compounds = []
    for item in inputs:
        data = json.loads(item.read_text())
        if data["schema_version"]!=3:
            raise ValueError(f"Expected schema-v3 preparation payload {item}")
        corrections.extend(data["corrections"])
        compounds.extend(data.get("compound_corrections", []))
    names = [node["name"] for node in corrections]
    if len(names)!=len(set(names)):
        raise ValueError("Duplicate preparation correction")
    required = {"tau_track_charged","tau_prepare","muon_selection_eta","muon_selection_vertex",
                "muon_selection_track_roles","muon_selection_after_calibration","electron_prepare","photon_prepare"}
    if not required.issubset(names):
        raise ValueError(f"Missing production preparation corrections: {sorted(required-set(names))}")
    combined = {"schema_version":3,"corrections":corrections}
    if compounds:
        combined["compound_corrections"] = compounds
    with output.open("x",encoding="utf-8") as handle:
        json.dump(combined,handle,allow_nan=False)
        handle.write("\n")


if __name__ == "__main__":
    typer.run(main)
