""" Functions for reading in PEAKS search results.
"""
import re

import pandas as pd

from ibench.constants import (
    ENGINE_SCORE_KEY,
    LABEL_KEY,
    PEPTIDE_KEY,
    SCAN_KEY,
    SOURCE_KEY,
)
from ibench.utils import remove_source_suffixes

# Define the relevant column names from PEAKS DB search results.
PEAKS_ACCESSION_KEY = 'Accession'
PEAKS_PEPTIDE_KEY = 'Peptide'
PEAKS_PTM_KEY = 'PTM'
PEAKS_SCAN_KEY = 'Scan'
PEAKS_SCORE_KEY = '-10lgP'
PEAKS_SOURCE_KEY = 'Source File'
PEAKS_RELEVANT_COLUMNS = [
    PEAKS_ACCESSION_KEY,
    PEAKS_PEPTIDE_KEY,
    PEAKS_PTM_KEY,
    PEAKS_SCAN_KEY,
    PEAKS_SCORE_KEY,
    PEAKS_SOURCE_KEY,
]


def read_single_peaks_data(df_loc, score_limit, hq_hits_only, filter_ptms):
    """ Function to read in PEAKS DB search results from a single file.

    Parameters
    ----------
    df_loc : str
        A location of PEAKS DB search results.

    Returns
    -------
    hits_df : pd.DataFrame
        A DataFrame of all search results properly formatted for Caravan.
    mods_dfs : pd.DataFrame
        A small DataFrame detailing the ptms found in the data.
    """
    peaks_df = pd.read_csv(df_loc, usecols=PEAKS_RELEVANT_COLUMNS)

    if filter_ptms:
        peaks_df = peaks_df[peaks_df[PEAKS_PEPTIDE_KEY].apply(lambda x : '(' not in x)]
    peaks_df[PEPTIDE_KEY] = peaks_df[PEAKS_PEPTIDE_KEY].apply(
        lambda x : re.sub(r'[^A-Za-z ]', '', x)
    )

    # Rename to match Caravan naming scheme.
    peaks_df = peaks_df.rename(columns={
        PEAKS_SCORE_KEY: ENGINE_SCORE_KEY,
    })

    if hq_hits_only:
        peaks_df = peaks_df[peaks_df[ENGINE_SCORE_KEY] > score_limit]

    peaks_df = peaks_df[peaks_df[PEAKS_ACCESSION_KEY].apply(lambda x : isinstance(x, str))]
    peaks_df[LABEL_KEY] = peaks_df[PEAKS_ACCESSION_KEY].apply(
        lambda x : -1 if isinstance(x, str) and '#DECOY#' in x else 1
    )

    if hq_hits_only:
        peaks_df = peaks_df[peaks_df[ENGINE_SCORE_KEY] > score_limit]
        peaks_df = peaks_df[peaks_df[LABEL_KEY] == 1]

    # Clean source and scan columns if required, add label.
    peaks_df[SOURCE_KEY] = peaks_df[PEAKS_SOURCE_KEY].apply(
        remove_source_suffixes
    )
    peaks_df[SCAN_KEY] = peaks_df[PEAKS_SCAN_KEY].apply(
        lambda x : x if isinstance(x, int) else int(x.split(':')[-1])
    )
    peaks_df = peaks_df[[
        SOURCE_KEY,
        SCAN_KEY,
        ENGINE_SCORE_KEY,
        PEPTIDE_KEY,
        LABEL_KEY,
    ]]


    return peaks_df
