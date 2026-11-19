"""Exercise the jet optional-branch projection through the real chunk packer."""

from contextlib import ExitStack
import unittest
from unittest.mock import patch

import awkward as ak
import numpy as np

from atlas_met_columnar import OWNERS, RAW_BRANCHES
from atlas_met_prepare import EVENT, JET, TRUTH, prepare_chunk
from met_demo import projected_fields

PS = JET + "PSFrac"
SAMPLES = JET + "EnergyPerSampling"
JET_OFFSETS = [0, 2, 2, 3]
SAMPLE_ROWS = [list(range(5)), list(range(6)), list(range(5))]


class JetProjection(unittest.TestCase):
    def project_and_pack(self, optional):
        numeric_terms = ak.Array([[0.], [0.], [0.]])
        raw = {
            **{EVENT + name: np.zeros(3) for name in (
                "runNumber", "eventTypeBitmask", "actualInteractionsPerCrossing",
                "averageInteractionsPerCrossing", "eventNumber", "mcChannelNumber",
            )},
            **{TRUTH + name: numeric_terms for name in ("mpx", "mpy", "sumet")},
            TRUTH + "source": ak.Array([[0x21000], [0x21000], [0x21000]]),
            TRUTH + "name": ak.Array([["NonInt"], ["NonInt"], ["NonInt"]]),
            JET + "NumTrkPt500": ak.Array([[[4], [2, 3]], [], [[1]]]),
            JET + "EMFrac": ak.Array([[0.5, 0.75], [], [0.6]]),
            JET + "Width": ak.Array([[0.13, 0.12], [], [0.11]]),
            **optional,
        }
        fields = projected_fields(raw)
        projected = {name: raw[name] for name in fields if name in raw}
        output = {short + ".o": np.zeros(4, dtype="i8") for short in OWNERS}
        output["jet.o"][:] = JET_OFFSETS
        with ExitStack() as stack:
            # Isolate the unrelated object preparation. Event/jet offsets,
            # optional masks and fallback packing use their production code.
            for name in ("primary_positions", "prepare_tau", "prepare_muon", "prepare_egamma"):
                stack.enter_context(patch("atlas_met_prepare." + name))
            prepare_chunk(output, projected, {}, None)
        return fields, output

    def test_present_psfrac_omits_the_fallback_even_if_stored(self):
        ps = ak.Array([[0.1, 0.3], [], [0.2]])
        fields, output = self.project_and_pack({PS: ps, SAMPLES: object()})
        self.assertEqual(fields[:len(RAW_BRANCHES)], RAW_BRANCHES)
        self.assertIn(PS, fields)
        self.assertNotIn(SAMPLES, fields)
        self.assertNotIn("jet.sample", output)
        self.assertNotIn("jet.sample.o", output)
        np.testing.assert_array_equal(output["jet.o"], JET_OFFSETS)
        np.testing.assert_array_equal(output["jet.ps.valid"], np.ones(3, dtype="u1"))
        np.testing.assert_array_equal(output["jet.ps.value"], np.array([0.1, 0.3, 0.2], dtype="f4"))

        _, without_stored_fallback = self.project_and_pack({PS: ps})
        self.assertEqual(output.keys(), without_stored_fallback.keys())
        for name in output:
            np.testing.assert_array_equal(output[name], without_stored_fallback[name])

    def test_missing_psfrac_preserves_jagged_sampling_fallback(self):
        samples = ak.Array([SAMPLE_ROWS[:2], [], SAMPLE_ROWS[2:]])
        fields, output = self.project_and_pack({SAMPLES: samples})
        self.assertNotIn(PS, fields)
        self.assertIn(SAMPLES, fields)
        np.testing.assert_array_equal(output["jet.o"], JET_OFFSETS)
        np.testing.assert_array_equal(output["jet.ps.valid"], np.zeros(3, dtype="u1"))
        np.testing.assert_array_equal(output["jet.ps.value"], np.zeros(3, dtype="f4"))
        np.testing.assert_array_equal(output["jet.sample.o"], [0, 5, 11, 16])
        np.testing.assert_array_equal(output["jet.sample"], np.concatenate(SAMPLE_ROWS).astype("f4"))

    def test_missing_both_keeps_the_original_missing_branch_failure(self):
        with self.assertRaises(KeyError) as error:
            self.project_and_pack({})
        self.assertEqual(error.exception.args, (SAMPLES,))
