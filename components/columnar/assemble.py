#!/usr/bin/env python3
"""Assemble a CPython 3.14 Linux wheel and its visible local uv script."""
from __future__ import annotations

import base64
import csv
from hashlib import sha256
from io import StringIO
import json
from pathlib import Path
import shutil
import subprocess
import sys
from zipfile import ZIP_DEFLATED, ZipFile

import typer

VERSION = "0.1.0"
TAG = "cp314-cp314-linux_x86_64"
DIST = "atlas_met_columnar"
PAYLOADS = ("maker.json", "nnjvt.json", "soft.json", "significance.json",
            "preparation.json", "EventFormat.json")
LIBRARIES = ("libatlas_met_columnar.so", "libatlas_met_maker_correction.so",
             "libatlas_met_maker_provider.so", "libmet_preparation_provider.so",
             "libcorrectionlib.so", "libstdc++.so.6", "libgcc_s.so.1")


def elf_deps(path):
    result = subprocess.run(["readelf","-d",str(path)],check=True,capture_output=True,text=True)
    deps = [line.split("[",1)[1].split("]",1)[0] for line in result.stdout.splitlines() if "(NEEDED)" in line]
    forbidden = ("libCore.so","libRIO.so","libHist.so","libTree.so","libMathCore","xAOD","Athena")
    if any(any(item in dep for item in forbidden) for dep in deps):
        raise ValueError(f"Nonportable native dependency in {path}: {deps}")
    if any("cvmfs" in line for line in result.stdout.splitlines() if "(RUNPATH)" in line or "(RPATH)" in line):
        raise ValueError(f"Install the native target with its relative RUNPATH first: {path}")
    return deps


def main(runtime: Path, output: Path):
    if sys.version_info[:2] != (3,14) or sys.platform != "linux":
        raise ValueError("Assemble from the selected CPython 3.14 Linux environment")
    source = Path(__file__).resolve().parent
    files = {
        "atlas_met_columnar.py": (source/"atlas_met_columnar.py").read_bytes(),
        "atlas_met_prepare.py": (source/"atlas_met_prepare.py").read_bytes(),
        "atlas_met_portable/__init__.py": (
            '"""Frozen ROOT-free MET provider and payloads."""\n'
            "from ctypes import CDLL, RTLD_GLOBAL\n"
            "from pathlib import Path\n"
            '_native = Path(__file__).parent / "native"\n'
            '_gcc = CDLL(str(_native / "libgcc_s.so.1"), mode=RTLD_GLOBAL)\n'
            '_cpp = CDLL(str(_native / "libstdc++.so.6"), mode=RTLD_GLOBAL)\n'
        ).encode(),
    }
    needed = {}
    for name in LIBRARIES:
        path = (runtime/"native"/name).resolve(strict=True)
        needed[name] = elf_deps(path)
        files["atlas_met_portable/native/"+name] = path.read_bytes()
    for name in PAYLOADS:
        files["atlas_met_portable/payloads/"+name] = (runtime/"payloads"/name).read_bytes()
    info = f"{DIST}-{VERSION}.dist-info"
    files[f"{info}/METADATA"] = (
        f"Metadata-Version: 2.3\nName: atlas-met-columnar\nVersion: {VERSION}\n"
        "Summary: Pinned AnalysisBase MET via CorrectionLib API v4\n"
        "Requires-Python: >=3.14,<3.15\nRequires-Dist: awkward==2.15.0\n"
        "Requires-Dist: uproot==5.7.7\nRequires-Dist: numpy\nRequires-Dist: hist\nRequires-Dist: typer\n"
    ).encode()
    files[f"{info}/WHEEL"] = f"Wheel-Version: 1.0\nGenerator: atlas-met-assemble\nRoot-Is-Purelib: false\nTag: {TAG}\n".encode()
    summary = {name:{"bytes":len(data),"sha256":sha256(data).hexdigest()} for name,data in files.items()}
    rows = StringIO(newline="")
    record = csv.writer(rows, lineterminator="\n")
    for name,data in sorted(files.items()):
        digest = base64.urlsafe_b64encode(sha256(data).digest()).rstrip(b"=").decode()
        record.writerow((name,"sha256="+digest,len(data)))
    record.writerow((f"{info}/RECORD","",""))
    files[f"{info}/RECORD"] = rows.getvalue().encode()
    output.mkdir()
    (output/"wheels").mkdir()
    wheel = output/"wheels"/f"{DIST}-{VERSION}-{TAG}.whl"
    with ZipFile(wheel,"w",ZIP_DEFLATED,compresslevel=6) as archive:
        for name,data in sorted(files.items()):
            archive.writestr(name,data)
    shutil.copyfile(source/"met_demo.py",output/"met_demo.py")
    with (output/"bundle.json").open("x") as handle:
        json.dump({"wheel":wheel.name, "wheel_bytes":wheel.stat().st_size,
                   "wheel_sha256":sha256(wheel.read_bytes()).hexdigest(),
                   "core_commit":"3dcfa20ec01fb40c4d5fec58383785519dec1c00",
                   "inputs":summary,"elf_needed":needed},handle,indent=2)
        handle.write("\n")
    print(wheel)


if __name__ == "__main__":
    typer.run(main)
