"""Load the frozen B bridge and run the off-production correlated diagnostic."""
from pathlib import Path
import ctypes
import json
import os
import sys

import typer

from run_corr_diagnostic import run

app = typer.Typer(add_completion=False)


@app.command()
def diagnostic(input_file: Path, yaml: Path, events: int) -> None:
    framework = Path(os.environ["ATLAS_MET_FRAMEWORK_PYTHON"]).resolve(strict=True).parent
    sys.path.insert(0, str(framework))
    from bridge_entry import audit_paths

    bridge = Path(os.environ["ATLAS_MET_FRAMEWORK_BRIDGE"]).resolve(strict=True)
    loaded = ctypes.CDLL(str(bridge), mode=ctypes.RTLD_GLOBAL)
    try:
        run(input_file, yaml, events, output_name="candidate")
    finally:
        with Path("binding-map.json").open("x") as output:
            json.dump(audit_paths(bridge, "full"), output, indent=2)
            output.write("\n")
    del loaded


if __name__ == "__main__":
    app()
