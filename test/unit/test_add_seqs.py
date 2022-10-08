""" Test suite for the iBench utils.
"""
import os
import random
import unittest

import numpy as np
import pandas as pd

from ibench.constants import (
    CANONICAL_KEY,
    CISSPLICED_KEY,
    TRANSPLICED_KEY,
)
from ibench.add_seqs import (
    add_canonical_seq,
    add_sequences,
    add_spliced_seq,
    get_insert_inds,
    remove_substring,
)
PROTEIN_SEQ = 'MITTAGESSENHIERSCHNELL'
PROTEIN_SEQ_2 = 'TESTINGSEQADD'
PROTEIN_SEQ_3 = 'YKRKLMNPQRS'
PROTEOME = [PROTEIN_SEQ, PROTEIN_SEQ_2, PROTEIN_SEQ_3]
PEPTIDE = 'HAPPY'

class TestAddSeqs(unittest.TestCase):
    def setUp(self):
        random.seed(42)
        np.random.seed(42)

    def test_add_canonical_seq(self):
        protein_seq, pep_idx = add_canonical_seq(PROTEIN_SEQ, PEPTIDE, [1])
        self.assertEqual(protein_seq, 'MHAPPYESSENHIERSCHNELL')
        self.assertEqual(pep_idx, 1)

        protein_seq, pep_idx = add_canonical_seq(PROTEIN_SEQ, PEPTIDE, [])
        self.assertEqual(protein_seq, 'MITTAGESSENHIERSCHNELLHAPPY')
        self.assertEqual(pep_idx, 22)

    def test_add_spliced_seq(self):
        protein_seq, insert_idx, splice_site_idx = add_spliced_seq(
            PROTEIN_SEQ, PEPTIDE, [2], splice_site_range=None
        )
        self.assertEqual(protein_seq, 'MIHLAPPYSSENHIERSCHNELL')
        self.assertEqual(insert_idx, 2)
        self.assertEqual(splice_site_idx, 1)

    def test_remove_substring(self):
        protein_seq = remove_substring(PROTEIN_SEQ, 'MITTAG')
        self.assertNotIn('MITTAG', protein_seq)

    def test_get_insert_inds(self):
        insert_inds = get_insert_inds([PROTEIN_SEQ, PROTEIN_SEQ, PROTEIN_SEQ], [PEPTIDE, PEPTIDE])
        self.assertEqual(insert_inds.tolist(), [2, 0])

    def test_add_sequences(self):
        proteome, tracking_df, modified_ids = add_sequences(
            PROTEOME,
            {CANONICAL_KEY: ['HAPPY'], CISSPLICED_KEY: ['TIRED'], TRANSPLICED_KEY: ['INRS']},
            None,
        )
        self.assertEqual(proteome, ['MITTAGETLIREDIERSCHNELLRS', 'TEININGSEQADD', 'YKRKLMNPQRHAPPY'])
        self.assertEqual(tracking_df['proteinIdx'].tolist(), [2, 0, -1])
        self.assertEqual(modified_ids, [2, 0, 1])
