import pandas as pd
import os
import unittest

from ibench.config import Config
from ibench.constants import (
    ENGINE_SCORE_KEY,
    PEPTIDE_KEY,
    LABEL_KEY,
    SCAN_KEY,
    SOURCE_KEY,
)

from ibench.input.maxquant import read_single_mq_data

EXPECTED_SINGLE_OUTPUT_COLUMNS = [
    SOURCE_KEY,
    SCAN_KEY,
    ENGINE_SCORE_KEY,
    PEPTIDE_KEY,
    LABEL_KEY,
]

class TestMaxQuant(unittest.TestCase):
    def setUp(self):
        cwd = os.getcwd()
        self.search_file_path = f'{cwd}/test/resources/max_quant_search.txt'

    def test_single_maxquant_read_all_filters(self):
        search_df = read_single_mq_data(
            self.search_file_path, 80, True, True
        )
        self.assertEqual(search_df.shape[0], 35)
        self.assertEqual(list(search_df.columns), EXPECTED_SINGLE_OUTPUT_COLUMNS)

    def test_single_maxquant_read_filter_keep_ptm(self):
        search_df = read_single_mq_data(
            self.search_file_path, 80, True, False
        )
        self.assertEqual(search_df.shape[0], 37)
        self.assertEqual(list(search_df.columns), EXPECTED_SINGLE_OUTPUT_COLUMNS)

    def test_single_maxquant_read_no_filter(self):
        search_df = read_single_mq_data(
            self.search_file_path, -1, False, False
        )
        self.assertEqual(search_df.shape[0], 408)
        self.assertEqual(list(search_df.columns), EXPECTED_SINGLE_OUTPUT_COLUMNS)

if __name__ == '__main__':
    unittest.main()


