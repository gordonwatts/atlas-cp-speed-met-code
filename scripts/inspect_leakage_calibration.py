"""Inventory the initialized egamma isolation leakage function expressions."""
from pathlib import Path
import json

import typer

app = typer.Typer(add_completion=False)
SOURCES = [
    "IsolationCorrections/v6/isolation_ptcorrections_rel22_mc20.root",
    "IsolationCorrections/v5/isolation_ptcorrections_rel20_2.root",
]
CLOSEBY = "ElectronIsolationSelection/v1/CaloDeltaRZ.root"


@app.command()
def inspect(output: Path) -> None:
    import ROOT

    ROOT.gSystem.Load("libPathResolverDict")
    report = {}
    for logical in SOURCES:
        resolved = str(ROOT.PathResolver.FindCalibFile(logical))
        root = ROOT.TFile.Open(resolved)
        forms = {}
        for key in root.GetListOfKeys():
            name = str(key.GetName())
            if not name.startswith("mean_f_topoetcone"):
                continue
            obj = key.ReadObj()
            identity = (
                str(obj.ClassName()), str(obj.GetTitle()), str(obj.GetExpFormula()),
                int(obj.GetNpar()), int(obj.GetNsave()) if hasattr(obj, "GetNsave") else None,
            )
            forms.setdefault(identity, []).append(name)
        report[logical] = {
            "resolved": resolved, "forms": [
                {
                    "class": cls, "title": title, "expression": expression,
                    "parameters": params, "saved_points": nsave, "keys": names,
                }
                for (cls, title, expression, params, nsave), names in forms.items()
            ],
        }
        root.Close()
    path = str(ROOT.PathResolver.FindCalibFile(CLOSEBY))
    root = ROOT.TFile.Open(path)
    obj = root.Get("hData")
    report[CLOSEBY] = {
        "resolved": path, "hData": {"class": str(obj.ClassName()), "dimensions": int(obj.GetDimension())},
    }
    root.Close()
    with output.open("x") as file:
        json.dump(report, file, indent=2, allow_nan=False)
        file.write("\n")


if __name__ == "__main__":
    app()
