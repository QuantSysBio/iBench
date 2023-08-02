""" Functions for reading in inSPIRE final assignments.
"""
import pandas as pd

from ibench.constants import (
    ENGINE_SCORE_KEY,
    LABEL_KEY,
    PEPTIDE_KEY,
    Q_VALUE_KEY,
    SCAN_KEY,
    SOURCE_KEY,
)

INSPIRE_SCORE_KEY = 'percolatorScore'
INSPIRE_Q_VALUE_KEY = 'qValue'

def read_single_inspire_data(
        psms_loc,
        q_value_limit,
        score_limit,
        hq_hits_only,
        label=1,
    ):
    """ Function to read an inSPIRE finalAssignments.csv file

    Parameters
    ----------
    psms_loc : str
        The location of the PSMs file.
    q_value_limit : float
        The cut off to be placed on the q-value.
    score_limit : float
        The cut off to be placed on the score.
    hq_hits_only : bool
        Indicator of whether we are selecting only high quality assignments.
    label : int
        Label of whether the PSMs are target (1) or decoy (-1).

    Returns
    -------
    inspire_df : pd.DataFrame
        The DataFrame of inSPIRE assignments.
    """
    inspire_df = pd.read_csv(psms_loc)

    if hq_hits_only:
        inspire_df = inspire_df[
            (inspire_df[INSPIRE_Q_VALUE_KEY] < q_value_limit) &
            (inspire_df[INSPIRE_SCORE_KEY] > score_limit)
        ]

    inspire_df = inspire_df.drop(ENGINE_SCORE_KEY, axis=1)
    inspire_df = inspire_df.rename(columns={
        INSPIRE_SCORE_KEY: ENGINE_SCORE_KEY,
        INSPIRE_Q_VALUE_KEY: Q_VALUE_KEY,
    })
    inspire_df[LABEL_KEY] = label
    inspire_df = inspire_df[[
        SOURCE_KEY,
        SCAN_KEY,
        ENGINE_SCORE_KEY,
        Q_VALUE_KEY,
        PEPTIDE_KEY,
        LABEL_KEY,
    ]]
    if inspire_df[PEPTIDE_KEY].iloc[0][1] == '.':
        inspire_df[PEPTIDE_KEY] = inspire_df[PEPTIDE_KEY].apply(
            lambda x : x[2:-2]
        )

    return inspire_df
