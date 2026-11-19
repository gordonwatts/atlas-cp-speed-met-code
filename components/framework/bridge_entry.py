#!/usr/bin/env python3
"""Load the released MET entry-point bridge inside the CPRun Python process."""
import ctypes
import json
import os
from pathlib import Path
import runpy


SYMBOLS = {
    "maker": "_ZN2CP11MetMakerAlg7executeERK12EventContext",
    "final": "_ZN2CP13MetBuilderAlg7executeERK12EventContext",
    "significance": "_ZN2CP18MetSignificanceAlg7executeERK12EventContext",
}


class DlInfo(ctypes.Structure):
    _fields_ = (
        ("filename", ctypes.c_char_p),
        ("base", ctypes.c_void_p),
        ("symbol", ctypes.c_char_p),
        ("address", ctypes.c_void_p),
    )


def audit_paths(bridge, stage):
    runtime = ctypes.CDLL(None)
    dladdr = runtime.dladdr
    dladdr.argtypes = (ctypes.c_void_p, ctypes.POINTER(DlInfo))
    dladdr.restype = ctypes.c_int
    operations = ("maker",) if stage == "maker_only" else tuple(SYMBOLS)
    records = {}
    for operation in operations:
        symbol = SYMBOLS[operation]
        address = ctypes.cast(getattr(runtime, symbol), ctypes.c_void_p)
        info = DlInfo()
        if dladdr(address, ctypes.byref(info)) != 1 or not info.filename:
            raise RuntimeError(f"Unable to verify live MET symbol {symbol}")
        owner = Path(os.fsdecode(info.filename)).resolve()
        if owner != bridge:
            raise RuntimeError(f"{symbol} resolved to {owner}, expected {bridge}")
        records[operation] = {"symbol": symbol, "owner": str(owner), "address": hex(address.value)}
    paths = {line.split()[-1] for line in Path("/proc/self/maps").read_text().splitlines()}
    maps = sorted(path for path in paths if path.startswith("/") and ".so" in path)
    return {
        "bridge": str(bridge), "stage": stage, "symbols": records, "mapped_libraries": maps,
        "worker_threads": len(list(Path("/proc/self/task").iterdir())),
        "cpu_affinity": sorted(os.sched_getaffinity(0)),
    }


def main():
    for option in ("LD_PRELOAD", "LD_DEBUG", "LD_DEBUG_OUTPUT"):
        os.environ.pop(option, None)
    bridge = Path(os.environ["ATLAS_MET_FRAMEWORK_BRIDGE"]).resolve(strict=True)
    stage = os.environ.get("ATLAS_MET_STAGE", "full")
    library = ctypes.CDLL(str(bridge), mode=ctypes.RTLD_GLOBAL)
    observer_path = os.environ.get("ATLAS_MET_FRAMEWORK_OBSERVER")
    observer = ctypes.CDLL(str(Path(observer_path).resolve(strict=True)), mode=ctypes.RTLD_GLOBAL) if observer_path else None
    runner = Path(os.environ["AnalysisBase_DIR"]) / "bin" / "CPRun.py"
    try:
        runpy.run_path(str(runner), run_name="__main__")
    finally:
        if os.environ.get("ATLAS_MET_DISPATCH_AUDIT"):
            evidence = audit_paths(bridge, stage)
            with Path("binding-map.json").open("x", encoding="utf-8") as handle:
                json.dump(evidence, handle, indent=2)
                handle.write("\n")


if __name__ == "__main__":
    main()
