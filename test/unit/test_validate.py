""" Test suite for the iBench validation of proteome.
"""
import os
import unittest

import numpy as np
import pandas as pd

from ibench.constants import (
    PEPTIDE_KEY,
)
from ibench.validate_assignments import check_assignment, validate_proteome

CAN_DF_ROW = {
    PEPTIDE_KEY: 'HAPPY',
    'proteinIdx': 2,
    'stratum': 'canonical',
    'frag1': None,
    'frag2': None,
}
CAN_DF_ROW_FAIL = {
    PEPTIDE_KEY: 'HARRY',
    'proteinIdx': 2,
    'stratum': 'canonical',
    'frag1': None,
    'frag2': None,
}
CIS_DF_ROW = {
    PEPTIDE_KEY: 'TIRED',
    'proteinIdx': 0,
    'stratum': 'cisspliced',
    'frag1': 'T',
    'frag2': 'IRED',
}
CIS_DF_ROW_FAIL = {
    PEPTIDE_KEY: 'HAPPY',
    'proteinIdx': 0,
    'stratum': 'cisspliced',
    'frag1': 'T',
    'frag2': 'IRED',
}
TRANS_DF_ROW = {
    PEPTIDE_KEY: 'INRS',
    'proteinIdx': -1,
    'stratum': 'transspliced',
    'frag1': 'IN',
    'frag2': 'IRED',
}
TRANS_DF_ROW_FAIL = {
    PEPTIDE_KEY: 'TIRED',
    'proteinIdx': -1,
    'stratum': 'transspliced',
    'frag1': 'T',
    'frag2': 'IRED',
}

MODIFIED_PROTEOME = ['MITTAGETLIREDIERSCHNELLRS', 'TEININGSEQADD', 'YKRKLMNPQRHAPPY']

class TestUtils(unittest.TestCase):
    def setUp(self):
        self.cwd = os.getcwd()
        self.hq_df = pd.read_csv(f'{self.cwd}/test/resources/high_confidence.csv')
        self.meta_df = pd.read_csv(f'{self.cwd}/test/resources/meta_df.csv')

    def test_validate_proteome_successful(self):
        n_invalid_entries = validate_proteome(self.hq_df, self.meta_df, 'test/resources/output', None, 25)
        self.assertEqual(n_invalid_entries, 0)

    def test_validate_proteome_failed(self):
        self.hq_df['stratum'] = pd.Series(['transspliced', 'canonical', 'cisspliced'])
        self.meta_df['stratum'] = pd.Series(['transspliced', 'canonical', 'cisspliced'])
        n_invalid_entries = validate_proteome(self.hq_df, self.meta_df, 'test/resources/output', None, 25)
        self.assertEqual(n_invalid_entries, 3)

    def test_check_assignment_can(self):
        assigned_successful = check_assignment(CAN_DF_ROW, MODIFIED_PROTEOME, True, None, 25)
        self.assertEqual(assigned_successful, True)
        assigned_fail = check_assignment(CAN_DF_ROW_FAIL, MODIFIED_PROTEOME, True, None, 25)
        self.assertEqual(assigned_fail, False)

    def test_check_assignment_cis(self):
        assigned_successful = check_assignment(CIS_DF_ROW, MODIFIED_PROTEOME, True, None, 25)
        self.assertEqual(assigned_successful, True)
        assigned_fail = check_assignment(CIS_DF_ROW_FAIL, MODIFIED_PROTEOME, True, None, 25)
        self.assertEqual(assigned_fail, False)

    def test_check_assignment_cis(self):
        assigned_successful = check_assignment(TRANS_DF_ROW, MODIFIED_PROTEOME, True, None, 25)
        self.assertEqual(assigned_successful, True)
        assigned_fail = check_assignment(TRANS_DF_ROW_FAIL, MODIFIED_PROTEOME, True, None, 25)
        self.assertEqual(assigned_fail, False)
