""" Function to extract confident hits.
"""
from copy import deepcopy
import difflib
import random

import pandas as pd
import numpy as np

from ibench.check_presence import generate_pairs
from ibench.constants import (
    CANONICAL_KEY,
    CISSPLICED_KEY,
    ENDC_TEXT,
    GT_SCAN_KEY,
    GT_SOURCE_KEY,
    IL_PEPTIDE_KEY,
    OKCYAN_TEXT,
    PEPTIDE_KEY,
    TRANSPLICED_KEY,
)
from ibench.input.mzml import process_mzml_files
from ibench.input.mgf import process_mgf_files
from ibench.input.search_results import generic_read_df

def _check_multiple_subseqs(peptide_list, max_group):
    """ Function to check if a peptide contains multiple subsequences,
        if so it shouldn't be assigned cisspliced.
    """
    for peptide1 in peptide_list:
        sub_pep = None
        for peptide2 in peptide_list:
            if peptide2 == peptide1:
                continue
            if peptide2 in peptide1:
                if sub_pep is None:
                    sub_pep = peptide2
                elif peptide2 in sub_pep:
                    sub_pep = peptide2
                elif sub_pep in peptide2:
                    continue
                else:
                    return random.choice([-1, max_group+1])
    return 0.0

def assign_strata(hq_df, config):
    """ Function to assign strata to each peptide in the ground truth dataset.

    Parameters
    ----------
    hq_df : pd.DataFrame
        A DataFrame of ground truth identified peptides.
    config : ibench.config.Config
        The Config object for the experiment.

    Returns
    -------
    hq_df : pd.DataFrame
        The input DataFrame with an additional column assigning accession strata.
    """
    print(
        OKCYAN_TEXT + f'\t{hq_df.shape[0]} High Quality Hits to embed in proteome.' + ENDC_TEXT
    )

    groups = hq_df['peptideGroup'].unique().tolist()
    shuffled_groups = list(np.random.choice(groups, size=len(groups), replace=False))
    hq_df['groupOrder'] = hq_df['peptideGroup'].apply(shuffled_groups.index)
    max_group = hq_df['groupOrder'].max()
    hq_df['forceCanonical'] = hq_df.groupby('groupOrder')[IL_PEPTIDE_KEY].transform(
        lambda x : _check_multiple_subseqs(list(x), max_group)
    )

    hq_df['groupOrder'] = hq_df[['groupOrder', 'forceCanonical']].apply(
        lambda x : x['forceCanonical'] if x['forceCanonical'] != 0 else x['groupOrder'],
        axis=1
    )
    hq_df = hq_df.sort_values(by='groupOrder')

    disc_cut_off = hq_df['groupOrder'].quantile(config.discoverable_fraction)
    ciss_cut_off = hq_df['groupOrder'].quantile(
        config.discoverable_fraction + config.cisspliced_fraction
    )

    if config.cisspliced_fraction > 0:
        groups = {
            CANONICAL_KEY: hq_df[
                hq_df['groupOrder'] < disc_cut_off
            ][IL_PEPTIDE_KEY].tolist(),
            CISSPLICED_KEY: hq_df[
                (hq_df['groupOrder'] >= disc_cut_off) &
                (hq_df['groupOrder'] < ciss_cut_off)
            ][IL_PEPTIDE_KEY].tolist(),
            TRANSPLICED_KEY: hq_df[hq_df['groupOrder'] >= ciss_cut_off][IL_PEPTIDE_KEY].tolist(),
        }
    else:
        groups = {
            CANONICAL_KEY: hq_df[
                hq_df['groupOrder'] < disc_cut_off
            ][IL_PEPTIDE_KEY].tolist(),
            CISSPLICED_KEY: [],
            TRANSPLICED_KEY: hq_df[hq_df['groupOrder'] >= disc_cut_off][IL_PEPTIDE_KEY].tolist(),
        }

    types_df = pd.DataFrame({
        IL_PEPTIDE_KEY: groups[CANONICAL_KEY] + groups[CISSPLICED_KEY] + groups[TRANSPLICED_KEY],
        'stratum': (
            [CANONICAL_KEY]*len(groups[CANONICAL_KEY]) +
            [CISSPLICED_KEY]*len(groups[CISSPLICED_KEY]) +
            [TRANSPLICED_KEY]*len(groups[TRANSPLICED_KEY])
        ),
    })

    hq_df = pd.merge(
        types_df,
        hq_df,
        how='inner',
        on=IL_PEPTIDE_KEY,
    )

    return hq_df

def check_closeness(peptide_1, peptide_2, similarity_cut_off):
    """ Function to check if two peptides are too similar to be assigned
        to different strata.
    """
    string_matcher = difflib.SequenceMatcher(
        None, peptide_1, peptide_2
    )
    _, __, size = string_matcher.find_longest_match(
        0, len(peptide_1), 0, len(peptide_2)
    )

    if len(peptide_1) <= len(peptide_2):
        if len(peptide_1) - size < similarity_cut_off:
            return True
        pairs = generate_pairs(peptide_1)
        for sr_1, sr_2 in pairs:
            if sr_1 in peptide_2 and sr_2 in peptide_2:
                return True

    if len(peptide_2) <= len(peptide_1):
        if len(peptide_2) - size < similarity_cut_off:
            return True
        pairs = generate_pairs(peptide_2)
        for sr_1, sr_2 in pairs:
            if sr_1 in peptide_1 and sr_2 in peptide_1:
                return True

    return False

def get_peptide_groups(peptides, similarity_cut_off):
    """ Function to group peptides into 3 separate strata so that peptides
        in different strata are not so different that it is impossible to
        embed them in the proteome.

    Parameters
    ----------
    peptides : list of str
        A list of the peptides to be added to the proteome.
    similarity_cut_off : int
        The minimum distance allowed for two peptides to be assigned to
        different strata.

    Returns
    -------
    peptide_groups : list of lists
        A list of groups of peptides which are so similar that they must
        be assigned to the same stratum.
    """
    peptide_groups = []
    for idx1, pep1 in enumerate(peptides):
        pep1_group = [pep1]
        for idx2, pep2 in enumerate(peptides):
            if idx1 == idx2:
                continue
            if check_closeness(pep1, pep2, similarity_cut_off):
                pep1_group.append(pep2)

        found_grp_inds = []
        for grp_idx, peptide_grp in enumerate(peptide_groups):
            for pep in pep1_group:
                if pep in peptide_grp:
                    found_grp_inds.append(grp_idx)

        if len(found_grp_inds) == 0:
            peptide_groups.append(pep1_group)
        elif len(found_grp_inds) == 1:
            peptide_groups[found_grp_inds[0]].extend(pep1_group)
        else:
            new_group = deepcopy(pep1_group)
            for found_id in sorted(list(set(found_grp_inds)), reverse=True):
                new_list = deepcopy(peptide_groups[found_id])
                new_group.extend(new_list)
                del peptide_groups[found_id]
            peptide_groups.append(new_group)

    return peptide_groups

def _get_group(peptide, peptide_groups):
    """ Helper function to get the peptide group of a peptide.
    """
    for idx, pep_grp in enumerate(peptide_groups):
        if peptide in pep_grp:
            return idx
    return -1

def extract_hq_hits(config):
    """ Function to extract the high quality PSMs from the search results provided.

    Parameters
    ----------
    config : ibench.config.Config
        The config object which controls the experiment.
    """
    # Read in all search results:
    target_df = generic_read_df(
        config,
        hq_hits_only=True,
    )

    # Create I/L redundant peptides
    target_df[IL_PEPTIDE_KEY] = target_df[PEPTIDE_KEY].apply(
        lambda x : x.replace('I', 'L')
    )
    target_df = target_df.drop_duplicates(subset=[IL_PEPTIDE_KEY])

    peptides = target_df[IL_PEPTIDE_KEY].tolist()
    peptide_groups = get_peptide_groups(peptides, config.closeness_cut_off)
    target_df['peptideGroup'] = target_df[IL_PEPTIDE_KEY].apply(
        lambda x : _get_group(x, peptide_groups)
    )

    hq_df = assign_strata(target_df, config)
    hq_df = hq_df.reset_index(drop=True)
    hq_df[GT_SOURCE_KEY] = f'ibenchGroundTruth_{config.identifier}'
    hq_df[GT_SCAN_KEY] = hq_df.index + 1

    if config.scan_format == 'mgf':
        hq_df = process_mgf_files(hq_df, config)
    else:
        hq_df = process_mzml_files(hq_df, config)

    hq_df = hq_df.sort_values(by=GT_SCAN_KEY)

    hq_df = hq_df.drop(
        ['peptideGroup', 'groupOrder', 'forceCanonical'],
        axis=1,
    )

    hq_df.to_csv(
        f'{config.output_folder}/high_confidence.csv',
        index=False,
    )
