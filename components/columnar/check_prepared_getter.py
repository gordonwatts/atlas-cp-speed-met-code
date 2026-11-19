#!/usr/bin/env python3
"""Check the registered v2 prepared momentum used by the production link adapter."""
from pathlib import Path

import numpy as np
import typer

from atlas_met_columnar import Native


def check(runtime: Path):
    native = runtime/"native"
    payloads = runtime/"payloads"
    bridge = Native(native/"libatlas_met_columnar.so")
    context = bridge.open(native/"libatlas_met_maker_provider.so",payloads/"maker.json",payloads/"nnjvt.json")
    try:
        bridge.configure(context,payloads/"soft.json",payloads/"significance.json",
                         native/"libmet_preparation_provider.so",payloads/"preparation.json")
        # ROOT PtEtaPhiM4D.h:150-189, etaMax.h:45 and AnalysisBase Photon_v1.cxx:46.
        # Native ParticleKind tags: jet=1, electron=2, photon=3, tau=4, muon=5.
        rows = {
            "pt": np.array([5.,5.,0.,0.,-5.,13.,0.]),
            "eta": np.array([0.,0.,22759.,-22759.,0.,0.,0.]),
            "phi": np.zeros(7),
            "mass": np.array([-3.,-6.,4.,4.,0.,0.,105.658367]),
            "particle_kind": np.array([1,1,1,4,3,3,5]),
        }
        got = bridge.prepare(context,"met_prepared_p4",rows)[:,2]
        expected = np.array([4.,0.,5.,5.,-5.,13.,105.658367])
        np.testing.assert_array_equal(got.view("u8"),expected.view("u8"))
    finally:
        bridge.lib.met_columnar_close(context)


if __name__ == "__main__":
    typer.run(check)
