""" Generic functions for reading in any search results.
"""
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import pandas as pd

from ibench.constants import (
    ENGINE_SCORE_KEY,
    HYDRO_INDEX_KEY,
    MASS_KEY,
    PEPTIDE_KEY,
    RESIDUE_WEIGHTS,
    SCAN_KEY,
    SEQ_LEN_KEY,
    SOURCE_KEY,
)
from ibench.input.mascot import read_single_mascot_data
from ibench.input.maxquant import read_single_mq_data
from ibench.input.peaks import read_single_peaks_data
from ibench.input.percolator import read_single_percolator_data

def _add_features(df_row):
    """ Function add features for to describe each PSM.

    Parameters
    ----------
    df_row : pd.Series
        A row of the PSM dataframe.

    Returns
    -------
    df_row : pd.Series
        The updated row containing calculated features.
    """
    peptide = df_row['peptide']
    df_row[HYDRO_INDEX_KEY] = ProteinAnalysis(peptide).gravy()
    df_row[SEQ_LEN_KEY] = len(peptide)
    df_row[MASS_KEY] = sum((RESIDUE_WEIGHTS[res] for res in peptide))

    return df_row

def combine_dataframes(all_search_results, hq_hits_only, config):
    """ Function to combine separate search results from different.

    Parameters
    ----------
    all_search_results : dict
        A dictionary containing identification groups mapped to lists of dataframes.
    hq_hits_only : bool
        Flag indicating whether
    config : ibench.config.Config
        The config object for the experiment.

    Returns
    -------
    combined_df : pd.DataFrame
        The identifications combined into a single DataFrame.
    """
    merged_result_list = []
    for id_group in all_search_results:
        if len(all_search_results[id_group]) > 1:
            id_df = all_search_results[id_group][0]
            for secondary_df in all_search_results[id_group][1:]:
                id_df = pd.merge(
                    id_df,
                    secondary_df[['source', SCAN_KEY, PEPTIDE_KEY]],
                    how='inner',
                    on=[SOURCE_KEY, SCAN_KEY, PEPTIDE_KEY],
                )
            merged_result_list.append(id_df)
        else:
            merged_result_list.append(all_search_results[id_group][0])
    combined_df = pd.concat(merged_result_list)
    combined_df = combined_df.sort_values(by=ENGINE_SCORE_KEY, ascending=False)

    if hq_hits_only:
        combined_df = combined_df.drop_duplicates(PEPTIDE_KEY)
        combined_df = combined_df[
            (combined_df['sequenceLength'] >= config.min_seq_len) &
            (combined_df['sequenceLength'] <= config.max_seq_len)
        ]

    return combined_df

def generic_read_df(config, hq_hits_only):
    """ Function to read in search results from any search engine.

    Parameters
    ----------
    config : ibench.config.Config
        The Config object for the experiment.
    hq_hits_only : bool
        Flag indicating whether to only accept PSMs above a given
        engine score/q-value threshold.

    Returns
    -------
    search_df : pd.DataFrame
        A DataFrame of search results.
    """
    all_results = {}
    for results in config.search_results:
        if results['searchEngine'] == 'mascot':
            search_df = read_single_mascot_data(
                results['resultsLocation'],
                results['qValueLimit'],
                results['scoreLimit'],
                hq_hits_only=hq_hits_only,
                filter_ptms=config.filter_ptms,
            )
        elif results['searchEngine'] == 'maxQuant':
            search_df = read_single_mq_data(
                results['resultsLocation'],
                results['scoreLimit'],
                hq_hits_only=hq_hits_only,
                filter_ptms=config.filter_ptms,
            )
        elif results['searchEngine'] == 'peaks':
            search_df = read_single_peaks_data(
                results['resultsLocation'],
                results['scoreLimit'],
                hq_hits_only=hq_hits_only,
                filter_ptms=config.filter_ptms,
            )
        elif results['searchEngine'] == 'percolator':
            search_df = read_single_percolator_data(
                results['resultsLocation'],
                results['qValueLimit'],
                results['scoreLimit'],
                hq_hits_only=hq_hits_only,
            )
        else:
            raise ValueError(f'Unknown Search Engine: {results["searchEngine"]}')

        search_df = search_df.apply(_add_features, axis=1)

        if results['identificationGroup'] in all_results:
            all_results[results['identificationGroup']].append(search_df)
        else:
            all_results[results['identificationGroup']] = [search_df]

    return combine_dataframes(all_results, hq_hits_only, config)
