""" Test suite for the iBench utils.
"""
import random
import unittest

import numpy as np


from ibench.check_presence import check_cis_present, find_cis_matched_splice_reactants, generate_pairs


PROTEIN_SEQ = 'MITTAGESSENHIERSCHNELL'
PEPTIDE = 'HAPPY'



class TestExtractHits(unittest.TestCase):
    def setUp(self):
        random.seed(42)
        np.random.seed(42)

    def test_check_cis(self):
        matched_reactants = find_cis_matched_splice_reactants(PROTEIN_SEQ, [('M', 'ILL'), ('MI', 'LL'), ('MIL', 'L')])
        self.assertEqual(matched_reactants, ['MI', 'LL'])
        matched_reactants = find_cis_matched_splice_reactants(PROTEIN_SEQ, [('M', 'ULL'), ('MU', 'LL'), ('MUL', 'L')])
        self.assertEqual(matched_reactants, None)

    def test_check_cis_present(self):
        cis_present = check_cis_present(PROTEIN_SEQ, 'MI', 'LL')
        self.assertEqual(cis_present, True)
        cis_present = check_cis_present(PROTEIN_SEQ, 'MU', 'LL')
        self.assertEqual(cis_present, False)

    def test_generate_pairs(self):
        splice_reactant_pairs = generate_pairs('MILL')
        self.assertEqual(splice_reactant_pairs, [('M', 'ILL'), ('MI', 'LL'), ('MIL', 'L')])
