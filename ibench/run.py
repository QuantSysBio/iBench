""" Main Script from which the whole program runs.
"""
from argparse import ArgumentParser
import random

import numpy as np

from ibench.config import Config
from ibench.constants import ENDC_TEXT, OKGREEN_TEXT
from ibench.download_data import download_data
from ibench.extract_hits import extract_hq_hits
from ibench.modify_db import modify_db
from ibench.performance import analyse_performance
from ibench.query_table import create_query_table

# pd.options.mode.chained_assignment = None

PIPELINE_OPTIONS = [
    'downloadExample',
    'createDB',
    'analysis',
]

random.seed(42)
np.random.seed(42)

def get_arguments():
    """ Function to collect command line arguments.

    Returns
    -------
    args : argparse.Namespace
        The parsed command line arguments.
    """
    parser = ArgumentParser(description='Benchmarking Mass Spectrometry Identification Methods.')

    parser.add_argument(
        '--config_file',
        help='Config file to be read from.'
    )

    parser.add_argument(
        '--pipeline',
        choices=PIPELINE_OPTIONS,
        required=True,
        help='What pipeline do you want to run?',
    )

    return parser.parse_args()


def main():
    """ Main function which executes the pipelines.
    """
    args = get_arguments()
    print(
        OKGREEN_TEXT +
        f'Running iBench {args.pipeline} pipeline!' +
        ENDC_TEXT
    )
    if args.pipeline != 'downloadExample':
        config = Config(args.config_file)
        config.validate(args.pipeline)

    if args.pipeline == 'downloadExample':
        download_data()

    # Creation of benchmarking datasets.
    if args.pipeline == 'createDB':
        extract_hq_hits(config)
        modify_db(config)

    # Analysis of results.
    if args.pipeline == 'analysis':
        create_query_table(config)
        analyse_performance(config)

if __name__ == '__main__':
    main()
