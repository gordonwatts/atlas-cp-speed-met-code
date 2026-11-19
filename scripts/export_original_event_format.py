"""Use xAOD dictionaries to witness link target metadata of the original file."""
from pathlib import Path
import json

import typer

app = typer.Typer(add_completion=False)


@app.command()
def export(input_file: Path, output_file: Path) -> None:
    import ROOT

    ROOT.xAOD.Init().ignore()
    if ROOT.gSystem.Load("libxAODEventFormatDict") < 0:
        raise RuntimeError("Could not load the release EventFormat dictionary.")
    source = ROOT.TFile.Open(str(input_file))
    meta = source.Get("MetaData")
    if meta.GetEntry(0) <= 0:
        raise RuntimeError("Could not read the original MetaData entry.")
    fmt = meta.EventFormatStreamDAOD_PHYSLITE
    fmt.toTransient()
    rows = {}
    for branch in source.Get("CollectionTree").GetListOfBranches():
        name = str(branch.GetName())
        element = fmt.get(name, True)
        if element:
            rows[name] = {
                "key": int(element.hash()),
                "class": str(element.className()),
                "parent": str(element.parentName()),
            }
    result = {
        "input_uuid": str(source.GetUUID().AsString()),
        "metadata_branch": "MetaData/EventFormatStreamDAOD_PHYSLITE",
        "top_level_eventformat_elements": rows,
    }
    assert rows["AnalysisJets"]["key"] == 452794649
    assert rows["AnalysisMuons"]["key"] == 980095599
    with output_file.open("x") as output:
        json.dump(result, output, indent=2, allow_nan=False)
        output.write("\n")
    source.Close()


if __name__ == "__main__":
    app()
