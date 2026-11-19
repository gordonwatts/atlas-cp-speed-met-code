#!/usr/bin/env python3
"""Run the M27 C++ fixture and check constructor calls in an untimed C trace."""
import ctypes as c
import json
from pathlib import Path

import typer


def main(diagnostic_library: Path, trace: Path | None = None,
         entries: int = typer.Option(16,min=0), chunk_size: int = typer.Option(7,min=1)):
    module = c.CDLL(str(diagnostic_library.resolve(strict=True)))
    module.met_columnar_error.argtypes = ()
    module.met_columnar_error.restype = c.c_char_p
    module.met_check_jet_ntrk_offsets.argtypes = ()
    module.met_check_jet_ntrk_offsets.restype = c.c_int
    if module.met_check_jet_ntrk_offsets() != 0:
        raise RuntimeError(module.met_columnar_error().decode())
    observed = []
    if trace is not None:
        with trace.open() as handle:
            for line in handle:
                entry,state,action,*values = line.rstrip("\n").split("\t")
                if action == "ntrkOffsetConstructorCalls":
                    observed.append((int(entry),int(state),*map(int,values)))
        expected = [(end-1,0,begin,end-begin,1)
                    for begin in range(0,entries,chunk_size) for end in [min(begin+chunk_size,entries)]]
        if observed != expected:
            raise ValueError(f"Unexpected jet.ntrk.o constructor calls {observed}, expected {expected}")
    print(json.dumps({"synthetic_cases":"passed","chunk_constructor_calls":observed}))


if __name__ == "__main__":
    typer.run(main)
