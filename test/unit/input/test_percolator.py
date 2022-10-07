import pandas as pd
import os
import unittest

from ibench.constants import (
    ENGINE_SCORE_KEY,
    PEPTIDE_KEY,
    LABEL_KEY,
    Q_VALUE_KEY,
    SCAN_KEY,
    SOURCE_KEY,
)

from ibench.input.percolator import read_single_percolator_data

EXPECTED_SINGLE_OUTPUT_COLUMNS = [
    SOURCE_KEY,
    SCAN_KEY,
    ENGINE_SCORE_KEY,
    Q_VALUE_KEY,
    PEPTIDE_KEY,
    LABEL_KEY,
]

class TeastPercolator(unittest.TestCase):
    def setUp(self):
        cwd = os.getcwd()
        self.search_file_path = f'{cwd}/test/resources/percolator.psms'

    def test_single_percolator_read_all_filters(self):
        search_df = read_single_percolator_data(
            self.search_file_path, 0.01, 0.0, True,
        )
        self.assertEqual(search_df.shape[0], 109)
        self.assertEqual(list(search_df.columns), EXPECTED_SINGLE_OUTPUT_COLUMNS)

    def test_single_percolator_read_no_filter(self):
        search_df = read_single_percolator_data(
            self.search_file_path, 1.01, -10.0, True,
        )
        self.assertEqual(search_df.shape[0], 165)
        self.assertEqual(list(search_df.columns), EXPECTED_SINGLE_OUTPUT_COLUMNS)

if __name__ == '__main__':
    unittest.main()
