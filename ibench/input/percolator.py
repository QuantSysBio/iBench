""" Methods for handling Percolator psms input format.
"""
import pandas as pd

from ibench.constants import (
    ENGINE_SCORE_KEY,
    SCAN_KEY,
    SOURCE_KEY,
)

def get_source_and_scan(df_row):
    """ Function to extract the source, scan, and charge from the PSMId.

    Parameters
    ----------
    df_row
    """
    psm_data = df_row['PSMId'].split('_')
    df_row[SOURCE_KEY] = '_'.join(psm_data[:-2])
    df_row[SCAN_KEY] = int(psm_data[-2])
    return df_row

def read_single_percolator_data(
        psms_loc,
        q_value_limit,
        score_limit,
        hq_hits_only,
    ):
    """ Function to read a Percolator .psms file

    Parameters
    ----------
    psms_loc : str
        The location of the PSMs file.
    q_value_limit : float
        The cut off we place on the q-value.
    """
    perc_df = pd.read_csv(psms_loc, sep='\t')

    if hq_hits_only:
        perc_df = perc_df[
            (perc_df['q-value'] < q_value_limit) &
            (perc_df['score'] > score_limit)
        ]

    perc_df = perc_df.apply(get_source_and_scan, axis=1)
    perc_df = perc_df.rename(columns={'score': ENGINE_SCORE_KEY})
    perc_df = perc_df[[
        SOURCE_KEY,
        SCAN_KEY,
        ENGINE_SCORE_KEY,
        'q-value',
        'peptide',
    ]]

    return perc_df
