""" Test suite for the iBench checking if peptide sequences are
    present as cisspliced peptides.
"""
import random
import unittest

import numpy as np

from ibench.check_presence import (
    check_cis_present,
    find_cis_matched_splice_reactants,
    generate_pairs,
)

PROTEIN_SEQ = 'MITTAGESSENHIERSCHNELL'

class TestCheckPresence(unittest.TestCase):
    def setUp(self):
        random.seed(42)
        np.random.seed(42)

    def test_generate_pairs(self):
        splice_reactant_pairs = generate_pairs('MILL')
        self.assertEqual(
            splice_reactant_pairs,
            [('M', 'ILL'), ('MI', 'LL'), ('MIL', 'L')],
        )

    def test_check_cis(self):
        matched_reactants = find_cis_matched_splice_reactants(
            PROTEIN_SEQ,
            [('M', 'ILL'), ('MI', 'LL'), ('MIL', 'L')]
        )
        self.assertEqual(matched_reactants, ['MI', 'LL'])

        matched_reactants = find_cis_matched_splice_reactants(
            PROTEIN_SEQ, [('M', 'ALL'), ('MA', 'LL'), ('MAL', 'L')]
        )
        self.assertEqual(matched_reactants, None)

    def test_check_cis_present(self):
        cis_present_mill = check_cis_present(PROTEIN_SEQ, 'MI', 'LL')
        self.assertEqual(cis_present_mill, True)

        cis_present_mall = check_cis_present(PROTEIN_SEQ, 'MA', 'LL')
        self.assertEqual(cis_present_mall, False)
