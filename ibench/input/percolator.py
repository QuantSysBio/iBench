""" Methods for handling Percolator psms input format.
"""
import pandas as pd

from ibench.constants import (
    ENGINE_SCORE_KEY,
    LABEL_KEY,
    PEPTIDE_KEY,
    SCAN_KEY,
    SOURCE_KEY,
    Q_VALUE_KEY
)
PERCOLATOR_SCORE_KEY = 'score'
PERCOLATOR_PSM_ID_KEY = 'PSMId'

def get_source_and_scan(df_row):
    """ Function to extract the source, scan, and charge from the PSMId.

    Parameters
    ----------
    df_row
    """
    psm_data = df_row[PERCOLATOR_PSM_ID_KEY].split('_')
    df_row[SOURCE_KEY] = '_'.join(psm_data[:-2])
    df_row[SCAN_KEY] = int(psm_data[-2])
    return df_row

def read_single_percolator_data(
        psms_loc,
        q_value_limit,
        score_limit,
        hq_hits_only,
        label=1,
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
            (perc_df[Q_VALUE_KEY] < q_value_limit) &
            (perc_df[PERCOLATOR_SCORE_KEY] > score_limit)
        ]

    perc_df = perc_df.apply(get_source_and_scan, axis=1)
    perc_df = perc_df.rename(columns={PERCOLATOR_SCORE_KEY: ENGINE_SCORE_KEY})
    perc_df[LABEL_KEY] = label
    perc_df = perc_df[[
        SOURCE_KEY,
        SCAN_KEY,
        ENGINE_SCORE_KEY,
        Q_VALUE_KEY,
        PEPTIDE_KEY,
        LABEL_KEY,
    ]]

    return perc_df
