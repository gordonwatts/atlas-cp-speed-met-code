"""Inventory the original CollectionTree's branch types and on-disk baskets."""
from pathlib import Path
import json

import typer
import uproot

app = typer.Typer(add_completion=False)


@app.command()
def inventory(input_file: Path, output_file: Path) -> None:
    with uproot.open(
        input_file, handler=uproot.source.file.MemmapSource,
        num_workers=1, num_fallback_workers=1,
    ) as original:
        tree = original["CollectionTree"]
        branches = {}
        for branch_name, branch in tree.items(recursive=True, full_paths=True):
            branches[branch_name] = {
                "type": branch.typename,
                "compressed_bytes": branch.compressed_bytes,
                "uncompressed_bytes": branch.uncompressed_bytes,
                "basket_count": branch.num_baskets,
            }
        report = {
            "input_uuid": str(original.file.uuid),
            "entry_count": tree.num_entries,
            "file_bytes": input_file.stat().st_size,
            "source": "uproot TBranch::fZipBytes/fTotBytes and basket inventory metadata",
            "branches": branches,
        }
    with output_file.open("x") as output:
        json.dump(report, output, indent=2)
        output.write("\n")


if __name__ == "__main__":
    app()
