import pandas as pd
import os
import unittest

from ibench.constants import (
    ENGINE_SCORE_KEY,
    LABEL_KEY,
    PEPTIDE_KEY,
    Q_VALUE_KEY,
    SCAN_KEY,
    SOURCE_KEY,
)

from ibench.input.mascot import read_single_mascot_data

EXPECTED_SINGLE_OUTPUT_COLUMNS = [
    SOURCE_KEY,
    SCAN_KEY,
    ENGINE_SCORE_KEY,
    Q_VALUE_KEY,
    PEPTIDE_KEY,
    LABEL_KEY,
]
class TestMascot(unittest.TestCase):
    def setUp(self):
        cwd = os.getcwd()
        self.target_file_path = f'{cwd}/test/resources/mascot_search.csv'
        self.decoy_file_path = f'{cwd}/test/resources/mascot_decoy_search.csv'

    def test_single_mascot_read_all_filters(self):
        search_df = read_single_mascot_data(
            self.target_file_path, 0.05, 10, True, True
        )
        self.assertEqual(search_df.shape[0], 29)
        self.assertEqual(list(search_df.columns), EXPECTED_SINGLE_OUTPUT_COLUMNS)

    def test_single_mascot_read_filter_keep_ptm(self):
        search_df = read_single_mascot_data(
            self.target_file_path, 0.05, 10, True, False
        )
        self.assertEqual(search_df.shape[0], 36)
        self.assertEqual(list(search_df.columns), EXPECTED_SINGLE_OUTPUT_COLUMNS)

    def test_single_mascot_read_no_filter(self):
        search_df = read_single_mascot_data(
            self.target_file_path, 1.0, -1, False, False
        )
        self.assertEqual(search_df.shape[0], 112)
        self.assertEqual(list(search_df.columns), EXPECTED_SINGLE_OUTPUT_COLUMNS)
    
    def test_single_mascot_read_decoy_filter(self):
        search_df = read_single_mascot_data(
            self.decoy_file_path, 1.0, -1, True, False
        )
        self.assertEqual(search_df.shape[0], 0)

    def test_single_mascot_read_decoy_no_filter(self):
        search_df = read_single_mascot_data(
            self.decoy_file_path, 1.0, -1, False, False
        )
        self.assertEqual(search_df.shape[0], 69)
        self.assertEqual(list(search_df.columns), EXPECTED_SINGLE_OUTPUT_COLUMNS)

if __name__ == '__main__':
    unittest.main()


