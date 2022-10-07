import os
import unittest

from ibench.config import Config
from ibench.constants import (
    ENGINE_SCORE_KEY,
    HYDRO_INDEX_KEY,
    LABEL_KEY,
    MASS_KEY,
    PEPTIDE_KEY,
    Q_VALUE_KEY,
    SCAN_KEY,
    SEQ_LEN_KEY,
    SOURCE_KEY,
)
from ibench.input.search_results import generic_read_df

EXPECTED_SINGLE_OUTPUT_COLUMNS_WITH_Q_VALUE = [
    SOURCE_KEY,
    SCAN_KEY,
    ENGINE_SCORE_KEY,
    Q_VALUE_KEY,
    PEPTIDE_KEY,
    LABEL_KEY,
    HYDRO_INDEX_KEY,
    SEQ_LEN_KEY,
    MASS_KEY,
]
EXPECTED_SINGLE_OUTPUT_COLUMNS_NO_Q_VALUE = [
    SOURCE_KEY,
    SCAN_KEY,
    ENGINE_SCORE_KEY,
    PEPTIDE_KEY,
    LABEL_KEY,
    HYDRO_INDEX_KEY,
    SEQ_LEN_KEY,
    MASS_KEY,
]
class TestSearchResults(unittest.TestCase):
    def setUp(self):
        self.cwd = os.getcwd()
        self.config = Config(f'{self.cwd}/test/resources/mock_config_file.yml')
        self.config.min_seq_len = 0
        self.config.max_seq_len = 100

    def test_generic_mascot_read(self):
        self.config.search_results = [
            {
                'resultsLocation': f'{self.cwd}/test/resources/mascot_search.csv',
                'identificationGroup': 0,
                'searchEngine': 'mascot',
                'qValueLimit': 0.05, 
                'scoreLimit': 10,
            },
        ]
        search_df = generic_read_df(self.config, True)

        self.assertEqual(search_df.shape[0], 28)
        self.assertEqual(search_df[PEPTIDE_KEY].nunique(), search_df.shape[0])
        self.assertEqual(list(search_df.columns), EXPECTED_SINGLE_OUTPUT_COLUMNS_WITH_Q_VALUE)

    def test_generic_peaks_read(self):
        self.config.search_results = [
            {
                'resultsLocation': f'{self.cwd}/test/resources/peaks_search.csv',
                'identificationGroup': 0,
                'searchEngine': 'peaks',
                'scoreLimit': 45,
            },
        ]
        search_df = generic_read_df(self.config, True)

        self.assertEqual(search_df.shape[0], 59)
        self.assertEqual(search_df[PEPTIDE_KEY].nunique(), search_df.shape[0])
        self.assertEqual(list(search_df.columns), EXPECTED_SINGLE_OUTPUT_COLUMNS_NO_Q_VALUE)

    def test_generic_mq_read(self):
        self.config.search_results = [
            {
                'resultsLocation': f'{self.cwd}/test/resources/max_quant_search.txt',
                'identificationGroup': 0,
                'searchEngine': 'maxQuant',
                'scoreLimit': 80,
            },
        ]
        search_df = generic_read_df(self.config, True)

        self.assertEqual(search_df.shape[0], 17)
        self.assertEqual(search_df[PEPTIDE_KEY].nunique(), search_df.shape[0])
        self.assertEqual(list(search_df.columns), EXPECTED_SINGLE_OUTPUT_COLUMNS_NO_Q_VALUE)

    def test_generic_percolator_read(self):
        self.config.search_results = [
            {
                'resultsLocation': f'{self.cwd}/test/resources/percolator.psms',
                'identificationGroup': 0,
                'searchEngine': 'percolator',
                'qValueLimit': 0.01, 
                'scoreLimit': 0,
            },
        ]
        search_df = generic_read_df(self.config, True)

        self.assertEqual(search_df.shape[0], 40)
        self.assertEqual(search_df[PEPTIDE_KEY].nunique(), search_df.shape[0])
        self.assertEqual(list(search_df.columns), EXPECTED_SINGLE_OUTPUT_COLUMNS_WITH_Q_VALUE)

if __name__ == '__main__':
    unittest.main()