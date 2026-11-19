"""Build and verify a transferable offline archive of the accepted MET wheel."""
from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
import shutil
from zipfile import ZIP_DEFLATED, ZipFile

from packaging.utils import parse_wheel_filename
import typer

app = typer.Typer(add_completion=False)

README = """# Standalone MET

Use CPython 3.14 and x86-64 Linux. The tested environment was AlmaLinux 9
under WSL2. The analysis executable includes its native CorrectionLib provider
and calibration payloads. Its input is the original PHYSLITE file identified
in input.txt. An exact local copy with the recorded ROOT UUID is also accepted.

From this extracted directory, run the five-state histogram with uv and the
included public dependency wheels. Use these thread limits to reproduce the
single-thread numerical environment.

```bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 BLIS_NUM_THREADS=1 RAYON_NUM_THREADS=1
uv run --offline --no-index --find-links public-wheels --python /path/to/python3.14 \\
    met_demo.py "$(cat input.txt)" full.npz
```

To create a separate nominal histogram, add `--nominal-only` and use a new
output filename. The separate correlated diagnostic uses
`--correlated-diagnostic`. The option `-n 2` reopens and recalculates the
original stream twice and checks that the final histogram repeats exactly.

The same local wheel can be installed with pip using the tested dependencies
and their exact hashes.

```bash
python3.14 -m venv .venv
.venv/bin/python -m pip install --no-index --find-links wheels \\
    --find-links public-wheels --require-hashes -r requirements-lock.txt
.venv/bin/python met_demo.py "$(cat input.txt)" full.npz
sha256sum -c files.sha256
```

Each Final histogram has unit weights, 100 bins from 0 to 500 GeV and both
flow bins. The NPZ contains labels, bin edges, counts, sum of weights and
sum of squared weights. Use a new output path so an existing result is
preserved. See bundle.json and identity.json for the compiled API identity,
payload hashes, original input and source revision.
"""


def digest(path: Path) -> str:
    with path.open("rb") as stream:
        return sha256(stream.read()).hexdigest()


def write_json(path: Path, content: dict) -> None:
    with path.open("x") as file:
        json.dump(content, file, indent=2)
        file.write("\n")


@app.command()
def assemble(portable: Path, source_head: Path, input_manifest: Path, output: Path) -> None:
    bundle = portable / "bundle"
    wheel, = sorted((bundle / "wheels").glob("*.whl"))
    manifest = json.loads((bundle / "bundle.json").read_text())
    if manifest["wheel_sha256"] != digest(wheel) or manifest["wheel_bytes"] != wheel.stat().st_size:
        raise ValueError("The observed portable wheel differs from its accepted manifest.")
    output.mkdir()
    stage = output / "bundle"
    stage.mkdir()
    (stage / "wheels").mkdir()
    shutil.copyfile(wheel, stage / "wheels" / wheel.name)
    for name in ("met_demo.py", "bundle.json", "source-layout.json"):
        shutil.copyfile(bundle / name, stage / name)
    shutil.copyfile(portable / "constraints.txt", stage / "constraints.txt")
    shutil.copyfile(input_manifest, stage / "input.txt")
    (stage / "README.md").write_text(README)
    public = stage / "public-wheels"
    public.mkdir()
    distribution, _, _, _ = parse_wheel_filename(wheel.name)
    locked = {}
    for dependency in sorted((portable / "wheelhouse").glob("*.whl")):
        name, version, _, _ = parse_wheel_filename(dependency.name)
        if name in locked:
            raise ValueError(f"Multiple tested wheel candidates for {name}")
        if name == distribution:
            if digest(dependency) != manifest["wheel_sha256"]:
                raise ValueError("The dependency wheelhouse contains a different MET wheel.")
        else:
            shutil.copyfile(dependency, public / dependency.name)
        locked[name] = f"{name}=={version} --hash=sha256:{digest(dependency)}"
    if distribution not in locked:
        raise ValueError("The tested wheelhouse does not include the local MET wheel.")
    (stage / "requirements-lock.txt").write_text("\n".join(locked[name] for name in sorted(locked)) + "\n")
    write_json(stage / "identity.json", {
        "source_commit": source_head.read_text().strip(),
        "script_sha256": digest(stage / "met_demo.py"),
        "wheel": wheel.name,
        "wheel_sha256": manifest["wheel_sha256"],
        "input_root_uuid": "79ee6652-311f-11f1-aa37-0445010abeef",
        "input_events": 20000,
        "python_abi": "cp314-cp314-linux_x86_64",
    })
    names = sorted(path.relative_to(stage) for path in stage.rglob("*") if path.is_file())
    hashes = {name.as_posix(): digest(stage / name) for name in names}
    with (stage / "files.sha256").open("x") as handle:
        for name, sha in hashes.items():
            handle.write(f"{sha}  {name}\n")
    archive = output / "met-standalone-cp314-linux-x86_64.zip"
    extracted = output / "verified-extraction"
    extracted.mkdir()
    with ZipFile(archive, "x", compression=ZIP_DEFLATED, compresslevel=6) as zipped:
        for name in [*names, Path("files.sha256")]:
            zipped.write(stage / name, name.as_posix())
    with ZipFile(archive) as zipped:
        invalid = zipped.testzip()
        if invalid:
            raise ValueError(f"Archive member failed its CRC: {invalid}")
        zipped.extractall(extracted)
    for name, sha in hashes.items():
        if digest(extracted / name) != sha:
            raise ValueError(f"Extracted archive member differs: {name}")
    write_json(output / "archive.json", {
        "archive": archive.name, "sha256": digest(archive), "bytes": archive.stat().st_size,
        "member_file_sha256": hashes, "all_extract_hashes_verified": True,
    })


if __name__ == "__main__":
    app()
