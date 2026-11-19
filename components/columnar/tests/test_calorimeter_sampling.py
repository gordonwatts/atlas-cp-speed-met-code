"""Compare packed cluster-sampling lookup with the original jagged path."""

from pathlib import Path
import os
import sys
import unittest
from unittest.mock import patch

import awkward as ak
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from atlas_met_prepare import CLUSTER, OWNERS, Data, Resolved, dense

SLOTS = (1, 2, 3, 5, 6, 7)
FIELDS = {"e_sampl": SLOTS, "eta_sampl": (1, 2, 5, 6)}


def branch(field):
    return CLUSTER + "AuxDyn." + field


def view(raw, rows, size=None, local=None):
    size = len(rows) if size is None else size
    local = np.arange(len(rows), dtype="i8") if local is None else np.asarray(local, dtype="i8")
    result = Resolved(raw, None, np.zeros(size, dtype="i8"), ())
    result.locations = [(CLUSTER, local, np.asarray(rows, dtype="i8"))]
    return result


def data(e=None, eta=None):
    raw = {branch(field): value for field, value in (("e_sampl", e), ("eta_sampl", eta)) if value is not None}
    output = {owner + ".o": np.zeros(1, dtype="i8") for owner in OWNERS}
    return Data(raw, output, {}, None, False)


def old_sampling(resolved, field, patterns, slot):
    output = np.zeros(resolved.size, dtype="f4")
    for target, local, rows in resolved.locations:
        source = ak.flatten(resolved.data.raw[target + "AuxDyn." + field], axis=1)[rows]
        relevant = (patterns[local] & (1 << slot)) != 0
        if not relevant.any():
            continue
        chosen = local[relevant]
        masks = patterns[chosen] & ((1 << slot) - 1)
        rank = np.zeros(len(chosen), dtype="i8")
        for bit in range(slot):
            rank += (masks >> bit) & 1
        selected = source[relevant]
        count = dense(ak.num(selected, axis=1), "i8")
        start = np.concatenate(([0], np.cumsum(count, dtype="i8")))
        good = rank < count
        output[chosen] = -999
        content = dense(ak.flatten(selected, axis=1), "f4")
        output[chosen[good]] = content[start[:-1][good] + rank[good]]
    return output


class SamplingProjection(unittest.TestCase):
    def assert_same(self, candidate, stock):
        self.assertEqual(candidate.dtype, np.dtype("f4"))
        np.testing.assert_array_equal(candidate.view("u4"), stock.view("u4"))

    def test_repeated_species_views_share_each_column_and_keep_rank_and_sentinels(self):
        original = data(
            ak.Array([[[100., 101., 105.], [205.], []], [], [[-1.5, 8.75, float("nan"), -0.0]]]),
            ak.Array([[[0.01, 0.04], [-1.0], []], [], [[0.02, 0.06, 0.08]]]),
        )
        electron = view(original, [0, 2, 3], size=4, local=[0, 2, 3])
        photon = view(original, [3, 1, 0, 3], size=5, local=[0, 1, 3, 4])
        e_pattern = np.array([0b10100111, 0xff, 0b00000010, 0b10100010], dtype="i8")
        p_pattern = np.array([0b10100010, 0b01000010, 0xff, 0b10100111, 0b00000010], dtype="i8")
        cases = [(linked, field, pattern, slot, old_sampling(linked, field, pattern, slot))
                 for linked, pattern in ((electron, e_pattern), (photon, p_pattern))
                 for field, slots in FIELDS.items() for slot in slots]
        with patch("atlas_met_prepare.ak.flatten", wraps=ak.flatten) as flatten:
            for linked, field, pattern, slot, stock in cases:
                with self.subTest(field=field, slot=slot, species=len(pattern)):
                    self.assert_same(linked.sampling(field, pattern, slot), stock)
        self.assertEqual(flatten.call_count, 2)
        self.assertEqual(set(original.sampling_columns), set(map(branch, FIELDS)))
        self.assertTrue(all(packed is not None for _, packed in original.sampling_columns.values()))
        np.testing.assert_array_equal(electron.sampling("e_sampl", e_pattern, 5),
                                      np.array([-999., 0., 0., 8.75], dtype="f4"))

    def test_noncontiguous_jagged_rows_and_source_offsets(self):
        primitive = ak.contents.NumpyArray(np.arange(10, dtype="f4"))
        for rows in (
            ak.contents.ListArray(ak.index.Index64(np.array([6, 0, 5, 2])),
                                  ak.index.Index64(np.array([8, 2, 5, 5])), primitive),
            ak.contents.ListOffsetArray(ak.index.Index64(np.array([2, 4, 4, 7, 8])), primitive),
        ):
            tree = ak.Array(ak.contents.ListOffsetArray(ak.index.Index64(np.array([0, 3, 4])), rows))
            original = data(tree)
            linked = view(original, [3, 0, 1, 0, 2])
            patterns = np.array([3, 34, 3, 36, 98], dtype="i8")
            for slot in SLOTS:
                self.assert_same(linked.sampling("e_sampl", patterns, slot),
                                 old_sampling(linked, "e_sampl", patterns, slot))

    def test_optional_and_empty_vectors_keep_original_behavior(self):
        for tree in (
            ak.Array([[None, [None, 2.]], [[4.]]]),
            ak.Array([[[], []], [[]]]),
        ):
            original = data(tree)
            linked = view(original, [0, 1, 2])
            patterns = np.array([2, 2, 2], dtype="i8")
            self.assert_same(linked.sampling("e_sampl", patterns, 1),
                             old_sampling(linked, "e_sampl", patterns, 1))
            self.assertIsNone(original.sampling_columns[branch("e_sampl")][1])
        unresolved = Resolved(data(), None, np.array([0, 0]), ())
        np.testing.assert_array_equal(unresolved.sampling("eta_sampl", np.array([0, 0]), 7),
                                      np.zeros(2, dtype="f4"))

    def test_missing_fields_and_out_of_range_gathers_fail_before_slot_filtering(self):
        for original, rows, error in (
            (data(), [0], KeyError),
            (data(), [], KeyError),
            (data(ak.Array([[[123.]]])), [1], IndexError),
        ):
            linked = view(original, rows)
            patterns = np.zeros(len(rows), dtype="i8")
            with self.assertRaises(error) as previous:
                old_sampling(linked, "e_sampl", patterns, 7)
            with self.assertRaises(error) as candidate:
                linked.sampling("e_sampl", patterns, 7)
            self.assertEqual(str(previous.exception), str(candidate.exception))

    @unittest.skipUnless(os.environ.get("ATLAS_MET_SAMPLING_PILOT"), "Set ATLAS_MET_SAMPLING_PILOT to original PHYSLITE")
    def test_original_numeric_layout_and_sample_values(self):
        import uproot

        with uproot.open(os.environ["ATLAS_MET_SAMPLING_PILOT"], handler=uproot.source.file.MemmapSource,
                         num_workers=1, num_fallback_workers=1) as file:
            self.assertEqual(str(file.file.uuid), "79ee6652-311f-11f1-aa37-0445010abeef")
            raw = file["CollectionTree"].arrays([branch(f) for f in (*FIELDS, "samplingPattern")],
                                               entry_stop=16, library="ak", how=dict)
        patterns = dense(ak.flatten(raw[branch("samplingPattern")], axis=1), "i8")
        self.assertGreater(len(patterns), 0)
        original = data(raw[branch("e_sampl")], raw[branch("eta_sampl")])
        linked = view(original, list(range(len(patterns))))
        for field, slots in FIELDS.items():
            for slot in slots:
                self.assert_same(linked.sampling(field, patterns, slot),
                                 old_sampling(linked, field, patterns, slot))
            self.assertIsNotNone(original.sampling_columns[branch(field)][1])
