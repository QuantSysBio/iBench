""" Functions for reading in MaxQuant search results.
"""
import pandas as pd

from ibench.constants import (
    ENGINE_SCORE_KEY,
    LABEL_KEY,
    PEPTIDE_KEY,
    SCAN_KEY,
    SOURCE_KEY,
)

# Define the relevant column names from MaxQuant search results.
MQ_DECOY_KEY = 'Reverse'
MQ_MODS_KEY = 'Modifications'
MQ_SCAN_KEY = 'Scan number'
MQ_SCORE_KEY = 'Score'
MQ_SEQ_KEY = 'Sequence'
MQ_SOURCE_KEY = 'Raw file'
MQ_RELEVANT_COLS = [
    MQ_DECOY_KEY,
    MQ_MODS_KEY,
    MQ_SCAN_KEY,
    MQ_SCORE_KEY,
    MQ_SEQ_KEY,
    MQ_SOURCE_KEY,
]

def read_single_mq_data(mq_data, score_limit, hq_hits_only, filter_ptms):
    """ Function to read in MaxQuant search results from a single file.

    Parameters
    ----------
    df_loc : str
        A location of MaxQuant search results.

    Returns
    -------
    hits_df : pd.DataFrame
        A DataFrame of all search results properly formatted for iBench.
    mods_dfs : pd.DataFrame
        A small DataFrame detailing the ptms found in the data.
    """
    mq_df = pd.read_csv(
        mq_data,
        sep='\t',
        usecols=MQ_RELEVANT_COLS,
    )

    if filter_ptms:
        mq_df = mq_df[mq_df[MQ_MODS_KEY] == 'Unmodified']

    # Rename to match iBench naming scheme.
    mq_df = mq_df.rename(
        columns={
            MQ_SCAN_KEY: SCAN_KEY,
            MQ_SCORE_KEY: ENGINE_SCORE_KEY,
            MQ_SEQ_KEY: PEPTIDE_KEY,
            MQ_SOURCE_KEY: SOURCE_KEY,
        }
    )

    mq_df[LABEL_KEY] = mq_df[MQ_DECOY_KEY].apply(
        lambda x : -1 if x == '+' else 1
    )

    mq_df = mq_df.dropna(subset=[PEPTIDE_KEY])

    if hq_hits_only:
        mq_df = mq_df[
            (mq_df[ENGINE_SCORE_KEY] > score_limit)
        ]
        mq_df = mq_df[mq_df[LABEL_KEY] == 1]

    mq_df = mq_df[[
        SOURCE_KEY,
        SCAN_KEY,
        ENGINE_SCORE_KEY,
        PEPTIDE_KEY,
        LABEL_KEY,
    ]]

    return mq_df
