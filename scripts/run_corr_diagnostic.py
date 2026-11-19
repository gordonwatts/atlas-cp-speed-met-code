"""Run the affecting-only correlated resolution state in a separate stock job."""
from __future__ import annotations

from pathlib import Path
import sys

import typer

app = typer.Typer(add_completion=False)


@app.command()
def run(input_file: Path, configuration: Path, events: int = typer.Option(16, min=1, max=100)) -> None:
    from AnalysisAlgorithmsConfig.EventLoopCPRunScript import EventLoopCPRunScript

    class CorrelatedMETDiagnostic(EventLoopCPRunScript):
        def modifyAlgSequence(self):
            (service,) = [alg for alg in self.algSeq if alg.type() == "CP::SystematicsSvc"]
            service.sigmaRecommended = 0
            service.systematicsList = ["", "MET_SoftTrk_ResoCorr"]

    sys.argv = [
        sys.argv[0], "--direct-driver", "--dump-full-config", "-i", str(input_file),
        "-t", str(configuration), "-o", "reference", "-e", str(events),
    ]
    CorrelatedMETDiagnostic().run()


if __name__ == "__main__":
    app()
