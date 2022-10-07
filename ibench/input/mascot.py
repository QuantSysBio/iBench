""" Functions for reading in Mascot search results.
"""
import pandas as pd

from ibench.constants import (
    ACCESSION_KEY,
    ENGINE_SCORE_KEY,
    LABEL_KEY,
    PEPTIDE_KEY,
    Q_VALUE_KEY,
    SCAN_KEY,
    SOURCE_KEY,
)

# Define separators within Mascot output and the names for relevant columns.
MASCOT_HEADER_MARKER = 'Header'
MASCOT_FILENAME_MARKER = 'Peak list data path'
MASCOT_HITS_START_MARKER = 'prot_hit_num'
MASCOT_NO_FIXED_MODS_MARKER = '""'
MASCOT_QUERIES_START_MARKER = 'Queries'
MASCOT_VAR_MODS_MARKER = 'Variable modifications'
MASCOT_FIXED_MODS_MARKER = 'Fixed modifications'
MASCOT_SEARCH_PARAMS_MARKER = 'Search Parameters'

MASCOT_MISS_KEY = 'pep_miss'
MASCOT_DECOY_KEY = 'prot_desc'
MASCOT_ENGINE_SCORE_KEY = 'pep_score'
MASCOT_PEPTIDE_KEY = 'pep_seq'
MASCOT_SCAN_TITLE_KEY = 'pep_scan_title'
MASCOT_PTM_SEQ_KEY = 'pep_var_mod_pos'
MASCOT_ACCESSION_KEY = 'prot_acc'
MASCOT_Q_VALUE_KEY = 'pep_expect'
REQUIRED_MASCOT_COLUMNS = [
    MASCOT_ACCESSION_KEY,
    MASCOT_DECOY_KEY,
    MASCOT_ENGINE_SCORE_KEY,
    MASCOT_PEPTIDE_KEY,
    MASCOT_PTM_SEQ_KEY,
    MASCOT_Q_VALUE_KEY,
    MASCOT_SCAN_TITLE_KEY,
]


def _get_mascot_file_metadata(csv_file):
    """ Function to extract metadata from a csv file of Mascot results.

    Parameters
    ----------
    csv_file : str
        The input address for a mascot search results file.

    Returns
    -------
    header_line : int
        The line number of the Header line.
    output_filename : str
        The name to be used for the output file of processing to ensure
        it matches with the processed Peaks file.
    hits_line : int
        The line number where the data on peptide matches starts.
    queries_line : int
        The line number where the query data starts (including Retention
        Time data).
    """
    with open(csv_file, 'r', encoding='UTF-8') as open_file:
        line_idx = 0
        line = open_file.readline()
        header_line = -1
        hits_line = -1
        queries_line = None
        search_params_line = -1

        while line:
            if header_line==-1 and MASCOT_HEADER_MARKER in line:
                header_line = line_idx
            if search_params_line == -1 and MASCOT_SEARCH_PARAMS_MARKER in line:
                search_params_line = line_idx
            if hits_line==-1 and MASCOT_HITS_START_MARKER in line:
                hits_line = line_idx
            if queries_line is None and MASCOT_QUERIES_START_MARKER in line:
                queries_line = line_idx
                break
            line = open_file.readline()
            line_idx += 1

    return header_line, hits_line, queries_line

def _skip_logic(idx, start_idx=0, end_idx=None):
    """ Simple function for identifying which lines to skip when reading a csv.

    Parameters
    ----------
    idx : int
        A line number in a file.
    start_idx : int
        The line number at which csv data starts.
    end_idx: int
        The line number at which csv data ends.

    Returns
    ------
    skip_line : bool
        A flag on whether to skip the line or not.
    """
    if end_idx is not None:
        if start_idx <= idx < end_idx:
            return False
    else:
        if idx >= start_idx:
            return False
    return True

def separate_scan_and_source(df_row):
    """ Function to separate source file and scan number (as well as retention time if
        Distiller format used) from mascots scan_title_format column.

    Parameters
    ----------
    df_row : pd.Series
        A row of the mascot search results DataFrame.
    scan_title_format : str or None
        A
    source_list : str or None
        A list of the source files if Mascot Distiller was used.

    Returns
    -------
    df_row : pd.Series
        The input row updated with new columns.
    """
    scan_title = df_row[MASCOT_SCAN_TITLE_KEY]
    df_row[SCAN_KEY] = int(scan_title.split('=')[-1].strip('~'))
    source = scan_title.split('File:')[-1]
    if source.startswith('~'):
        source = source.split('~')[1]
    else:
        source = source.split(', ')[0]
    if source.endswith('.raw'):
        source = source[:-4]
    df_row[SOURCE_KEY] = source

    return df_row

def _read_mascot_dfs(
        csv_filename,
        hits_line,
        queries_line,
        q_value_limit,
        score_limit,
        hq_hits_only,
        filter_ptms,
    ):
    """ Function to read the various sections of a Mascot output file.

    Parameters
    ----------
    csv_filename : str
        The path to the file containing all of the Mascot output.
    hits_line : int
        The line where the hits (PSMs) begin.
    queries_line : int
        The line where the query data begins.
    q_value_limit : float
        The q-value cut off above which PSMs will be filtered.
    score_limit : float
        The Mascot ion score cut off below which PSMs will be filtered.
    hq_hits_only : bool
        Boolean flag on whehther to apply score and q-value filters and drop decoys.
    filter_ptms : bool
        Boolean flag on whether to remove PTM modified peptides.

    Returns
    -------
    hits_df : pd.DataFrame
        A DataFrame of PSMs.
    """
    hits_df = pd.read_csv(
        csv_filename,
        skiprows=lambda idx : _skip_logic(idx, hits_line, queries_line),
        usecols=REQUIRED_MASCOT_COLUMNS
    )
    if filter_ptms:
        hits_df = hits_df[hits_df[MASCOT_PTM_SEQ_KEY].isna()]

    # Rename to match Caravan naming scheme.
    hits_df = hits_df.rename(
        columns={
            MASCOT_ACCESSION_KEY: ACCESSION_KEY,
            MASCOT_ENGINE_SCORE_KEY: ENGINE_SCORE_KEY,
            MASCOT_Q_VALUE_KEY: Q_VALUE_KEY,
            MASCOT_PEPTIDE_KEY: PEPTIDE_KEY,
            MASCOT_MISS_KEY: 'missedCleavages',
        }
    )
    hits_df = hits_df.sort_values(by=ENGINE_SCORE_KEY, ascending=False)
    hits_df = hits_df.drop_duplicates(MASCOT_SCAN_TITLE_KEY)

    hits_df[LABEL_KEY] = hits_df[MASCOT_DECOY_KEY].apply(
        lambda x : -1 if 'Reversed' in x or 'Random' in x else 1
    )
    if hq_hits_only:
        hits_df = hits_df[
            (hits_df[Q_VALUE_KEY] < q_value_limit) &
            (hits_df[ENGINE_SCORE_KEY] > score_limit)
        ]
        hits_df = hits_df[hits_df[LABEL_KEY] == 1]

    # Clean source and scan columns if required, add label.
    hits_df = hits_df.apply(separate_scan_and_source, axis=1)

    return hits_df

def read_single_mascot_data(input_filename, q_value_limit, score_limit, hq_hits_only, filter_ptms):
    """ Function to read in mascot search results from a single file.

    Parameters
    ----------
    input_filename : str
        A location of mascot search results.

    Returns
    -------
    hits_df : pd.DataFrame
        A DataFrame of all search results properly formatted for Caravan.
    mods_dfs : pd.DataFrame
        A small DataFrame detailing the ptms found in the data.
    """
    # Original File contains both peptide hits and queries, search file to find separators.
    header_line, hits_line, queries_line = _get_mascot_file_metadata(
        input_filename
    )

    if header_line == -1:
        raise ValueError('No header line found in Mascot search results.')

    # Read separate files as pandas dataframes.
    hits_df = _read_mascot_dfs(
        input_filename,
        hits_line,
        queries_line,
        q_value_limit,
        score_limit,
        hq_hits_only,
        filter_ptms,
    )

    if hits_df.shape[0]:
        hits_df = hits_df[[
            'source',
            'scan',
            ENGINE_SCORE_KEY,
            Q_VALUE_KEY,
            'peptide',
            LABEL_KEY,
        ]]

    return hits_df
