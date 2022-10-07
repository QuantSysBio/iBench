""" Test suite for the iBench utils.
"""
import os
import unittest

import numpy as np
import pandas as pd

from ibench.constants import (
    CANONICAL_KEY,
    CHARGE_KEY,
    CISSPLICED_KEY,
    INTENSITIES_KEY,
    ION_OFFSET,
    MZS_KEY,
    PEPTIDE_KEY,
    PROTON,
    TRANSPLICED_KEY,
)
from ibench.utils import calculate_ms2_feats, compute_potential_mzs, get_matches, get_pepitde_strata, remove_source_suffixes

EXPECTED_MZS = [71.037114, 174.046299, 289.073242, 418.115835]
EXPECTED_REVERSED_MZS = [147.068414, 276.111007, 391.13795, 494.147135]
GROUP_MASSES_B = np.array(EXPECTED_MZS) + ION_OFFSET['b']
TEST_DF_ROW = {
    PEPTIDE_KEY: 'ACDEF',
    CHARGE_KEY: 1,
    MZS_KEY: np.array([
        ION_OFFSET['b'] + 71.04 + PROTON,
        ION_OFFSET['y'] + 276.1113 + PROTON,
        ION_OFFSET['a'] + 333.22 + PROTON,
    ]),
    INTENSITIES_KEY: np.array([
        200_000,
        500_000,
        300_000,
    ]),
}

class TestUtils(unittest.TestCase):
    def setUp(self):
        self.cwd = os.getcwd()
        self.hq_df = pd.read_csv(f'{self.cwd}/test/resources/high_confidence.csv')

    def test_remove_source_suffixes(self):
        mgf_removed = remove_source_suffixes('test.mgf')
        self.assertEqual(mgf_removed, 'test')
        mzml_removed = remove_source_suffixes('test.mzML')
        self.assertEqual(mzml_removed, 'test')
        raw_removed = remove_source_suffixes('test.raw')
        self.assertEqual(raw_removed, 'test')

    def test_get_matches(self):
        matched_location, matched_index = get_matches(TEST_DF_ROW, GROUP_MASSES_B, 1, 0.02)
        self.assertEqual(matched_location, [1])
        self.assertEqual(matched_index, [0])

    def test_get_pepitde_strata(self):
        peptide_strata = get_pepitde_strata(self.hq_df)
        self.assertEqual(len(peptide_strata[CANONICAL_KEY]), 1)
        self.assertEqual(len(peptide_strata[CISSPLICED_KEY]), 1)
        self.assertEqual(len(peptide_strata[TRANSPLICED_KEY]), 1)

    def test_calculate_ms2_feats(self):
        updated_df_row = calculate_ms2_feats(TEST_DF_ROW, 0.02)
        self.assertEqual(updated_df_row['ms2Coverage'], 0.5)
        self.assertEqual(updated_df_row['signalToNoise'], 0.7)

    def test_compute_potential_mzs(self):
        potential_mzs = compute_potential_mzs('ACDEF', reverse=False)
        rev_potential_mzs = compute_potential_mzs('ACDEF', reverse=True)
        self.assertEqual(potential_mzs.tolist(), EXPECTED_MZS)
        self.assertEqual(rev_potential_mzs.tolist(), EXPECTED_REVERSED_MZS)
