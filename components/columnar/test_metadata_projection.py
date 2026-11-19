"""Verify exact-name projection on direct, nested and split ROOT metadata."""

import unittest

from atlas_met_columnar import RAW_BRANCHES
from atlas_met_prepare import ADDITIONAL, JET
from met_demo import projected_fields, read_branch_names

PS = JET + "PSFrac"
SAMPLES = JET + "EnergyPerSampling"
TRUTH_LINK = "AnalysisTauJetsAuxDyn.truthParticleLink"
ABSENT = "AnalysisTauJetsAuxDyn.TESCompatibility"


class BranchMetadata:
    def __init__(self, paths):
        self.paths = paths
        self.key_calls = []
        self.contains_calls = []

    def keys(self, *, recursive=True, full_paths=True):
        self.key_calls.append((recursive, full_paths))
        names = (path for path in self.paths if recursive or len(path) == 1)
        return ["/".join(path) if full_paths else path[-1] for path in names]

    def __contains__(self, name):
        self.contains_calls.append(name)
        return any(path[-1] == name for path in self.paths)


class MetadataProjection(unittest.TestCase):
    def compare_with_direct_lookup(self, paths):
        metadata = BranchMetadata(paths)
        previous = projected_fields(metadata)
        previous_absent = [name for name in ADDITIONAL if name not in metadata]
        previous_calls = metadata.contains_calls.copy()

        available = read_branch_names(metadata)
        fields = projected_fields(available)
        absent = [name for name in ADDITIONAL if name not in available]
        self.assertEqual(metadata.key_calls, [(True, False)])
        self.assertEqual(metadata.contains_calls, previous_calls)
        self.assertEqual(fields, previous)
        self.assertEqual(absent, previous_absent)
        self.assertEqual(fields[:len(RAW_BRANCHES)], RAW_BRANCHES)
        return available, fields, absent

    def test_recursive_exact_dotted_names_and_input_order(self):
        # Metadata order differs from the projection order. A split parent
        # must itself be selected, and its children remain distinct names.
        available, fields, absent = self.compare_with_direct_lookup([
            (SAMPLES,), ("EventInfoAuxDyn.runNumber",),
            ("SomeParent.", PS), ("AnotherParent.", PS),
            (TRUTH_LINK,), (TRUTH_LINK, TRUTH_LINK + ".m_persKey"),
            (TRUTH_LINK, TRUTH_LINK + ".m_persIndex"),
            (ABSENT + "_other",), ("MET_TruthAuxDyn.mpx",),
        ])
        self.assertIn(PS, fields)
        self.assertNotIn(SAMPLES, fields)
        self.assertNotIn(SAMPLES, absent)
        self.assertIn(TRUTH_LINK, fields)
        self.assertIn(TRUTH_LINK + ".m_persIndex", available)
        self.assertNotIn(TRUTH_LINK + ".m_persIndex", fields)
        self.assertIn(ABSENT, absent)
        self.assertEqual(fields[len(RAW_BRANCHES):], [
            "EventInfoAuxDyn.runNumber", "MET_TruthAuxDyn.mpx", PS, TRUTH_LINK,
        ])

    def test_nested_fallback_retained_when_psfrac_absent(self):
        _, fields, absent = self.compare_with_direct_lookup([
            ("SomeParent.", SAMPLES), (PS + "_other",),
        ])
        self.assertIn(SAMPLES, fields)
        self.assertIn(PS, absent)
