"""Freeze minimal stock production and stored-core control YAML from the observer."""
from pathlib import Path

import typer
import yaml

app = typer.Typer(add_completion=False)


@app.command()
def prepare(observer_yaml: Path, production_yaml: Path, stored_yaml: Path) -> None:
    configuration = yaml.safe_load(observer_yaml.read_text())
    configuration.pop("PrintConfiguration")
    configuration["Output"] = {
        "metVars": ["AnaMET_%SYS%.met -> finalMet_met_%SYS% type=float"],
    }
    with production_yaml.open("x") as output:
        yaml.safe_dump(configuration, output, sort_keys=False)
    stored = {
        "CommonServices": {"runSystematics": False, "systematicsHistogram": "systematicsList"},
        "Output": {
            "vars": [
                f"MET_Core_AnalysisMET.{component} -> core_{component} type=float metTerm=PVSoftTrkCore"
                for component in ("mpx", "mpy")
            ]
        },
    }
    with stored_yaml.open("x") as output:
        yaml.safe_dump(stored, output, sort_keys=False)


if __name__ == "__main__":
    app()
