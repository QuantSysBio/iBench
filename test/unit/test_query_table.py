""" Test suite for the iBench query table functions.
"""
import unittest

from ibench.constants import (
    CANONICAL_KEY,
    CISSPLICED_KEY,
)
from ibench.query_table import remap_to_proteome

PROTEOME = ['MITTAGETLIREDIERSCHNELLRS', 'TEININGSEQADD', 'YKRKLMNPQRHAPPY']


class TestQueryTable(unittest.TestCase):
    def test_remap_to_proteome_can(self):
        stratum = remap_to_proteome('HAPPY', PROTEOME, 25)
        self.assertEqual(stratum, CANONICAL_KEY)

    def test_remap_to_proteome_cis(self):
        stratum = remap_to_proteome('TIRED', PROTEOME, 25)
        self.assertEqual(stratum, CISSPLICED_KEY)

    def test_remap_to_proteome_trans(self):
        stratum = remap_to_proteome('INRS', PROTEOME, 25)
        self.assertEqual(stratum, 'unknown')

