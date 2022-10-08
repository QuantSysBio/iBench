""" Functions for adding sequences to the modified proteome.
"""
import random

import numpy as np
import pandas as pd

from ibench.constants import AMINO_ACIDS, CANONICAL_KEY, CISSPLICED_KEY, TRANSPLICED_KEY

def add_canonical_seq(protein_seq, peptide, choices):
    """ Function to add a canonical peptide to a protein at a random location.

    Parameters
    ----------
    protein_seq : str
        A full protein sequence.
    peptide : str
        The peptide sequence to be embedded in the protein.
    choices : list of int
        A list of the the possible indices at which the peptide can add.

    Returns
    -------
    protein_seq : str
        The modified protein sequence.
    position : int
        The index at which the peptide was added.
    """
    if len(choices):
        position = random.choice(choices)
        protein_seq = (
            protein_seq[:position] +
            peptide +
            protein_seq[position + len(peptide):]
        )
    else:
        position = len(protein_seq)
        protein_seq += peptide
    return protein_seq, position

def add_spliced_seq(protein_seq, peptide, choices, splice_site_range=None):
    """ Function to add a spliced peptide to a protein sequence.

    Parameters
    ----------
    protein_seq : str
        The protein sequence.
    peptide : str
        The peptide sequence.
    choices : list of int
        The possible indices where we can edit the protein.
    splice_site_range : list of int or None
        The indices at which the peptide can be split.

    Returns
    -------
    protein_seq : str
        The modified protein sequence.
    position : int
        The index at which the first peptide fragment was added.
    splice_site : int
        The index on the peptide was split.
    """
    pep_len = len(peptide)
    if splice_site_range is None:
        splice_site = random.randint(1, pep_len - 1)
    else:
        splice_site = random.choice(splice_site_range)
    int_seq = ''.join(
        [random.choice(AMINO_ACIDS) for _ in range(random.randint(1, 26))]
    )

    insert_string = peptide[:splice_site] + int_seq + peptide[splice_site:]

    if len(choices):
        position = random.choice(choices)
        protein_seq = (
            protein_seq[:position] +
            insert_string +
            protein_seq[position+pep_len:]
        )
    else:
        position = len(protein_seq)
        protein_seq += (
            peptide[:splice_site] +
            int_seq +
            peptide[splice_site:]
        )

    return protein_seq, position, splice_site

def remove_substring(protein, sub_string):
    """ Helper function to remove a substring from a protein.

    Parameters
    ----------
    protein : str
        The protein from which we are removing an entry.
    sub_string : str
        The substring being removed.

    Returns
    -------
    protein : str
        The input protein with the sub_string removed.
    """
    replacement = ''.join(
        [random.choice(AMINO_ACIDS) for _ in range(len(sub_string))]
    )
    return protein.replace(sub_string, replacement)


def get_insert_inds(proteome, peptide_list):
    """ Function to randomly chose a protein for peptide addition.

    Parameters
    ----------
    proteome : list of str
        A list of the proteins found in the proteome.
    peptide_list : list of str
        A list of the peptides to be added.

    Returns
    -------
    choices : list of int
        A list of the indices where peptides will be added.
    """
    n_proteins = len(proteome)
    return np.random.choice(
        range(n_proteins),
        size=len(peptide_list),
    )

def get_trans_insert_inds(proteome, peptide_list):
    """ Function to randomly chose two proteins for peptide fragments addition.

    Parameters
    ----------
    proteome : list of str
        A list of the proteins found in the proteome.
    peptide_list : list of str
        A list of the peptides to be added.

    Returns
    -------
    insert_inds_1 : list of int
        A list of the indices where the first fragment will be added.
    insert_inds_2 : list of int
        A list of the indices where the second fragment will be added.
    """
    n_proteins = len(proteome)
    insert_inds_1 = []
    insert_inds_2 = []
    for _ in peptide_list:
        trans_frag_inds = np.random.choice(
            range(n_proteins),
            size=2,
            replace=False,
        )
        insert_inds_1.append(trans_frag_inds[0])
        insert_inds_2.append(trans_frag_inds[1])

    return insert_inds_1, insert_inds_2


def add_sequences(proteome, peptide_strata, enzyme):
    """ Function to add sequences to the proteome.

    Parameters
    ----------
    proteome : list of str
        The list protein sequences to be modified.
    peptide_strata : dict
        The peptides in the ground truth data.
    enzyme : str
        Trypsin or unspecific.

    Returns
    -------
    modified_proteome : list of str
        The input proteome modified to include the ground truth peptides.
    peptide_df : pd.DataFrame
        A DataFrame containing details of where peptides were added.
    modified_ids : list of str
        A list of the indices which were modified.
    """
    modified_ids = {}
    peptide_data = []
    insert_inds = get_insert_inds(proteome, peptide_strata[CANONICAL_KEY])
    for peptide, prot_idx in zip(peptide_strata[CANONICAL_KEY], insert_inds):
        if prot_idx in modified_ids:
            choices = []
            for i in range(len(proteome[prot_idx])):
                idx_ok = True
                for mod_pos in modified_ids[prot_idx]:
                    if abs(i-mod_pos) < 20:
                        idx_ok = False
                if idx_ok:
                    choices.append(i)
        else:
            modified_ids[prot_idx] = []
            choices = list(range(len(proteome[prot_idx])))

        if enzyme is None:
            proteome[prot_idx], mod_pos = add_canonical_seq(
                proteome[prot_idx], peptide, choices
            )
        elif enzyme == 'trypsin':
            proteome[prot_idx], mod_pos = add_canonical_seq(
                proteome[prot_idx], 'K' + peptide, choices
            )
        else:
            raise ValueError('Unrecognised enzyme. Allowed values are None, trypsin.')
        modified_ids[prot_idx].append(mod_pos)
        peptide_data.append({
            'peptide': peptide,
            'proteinIdx': prot_idx,
            'stratum': CANONICAL_KEY,
            'frag1': None,
            'frag2': None,
        })

    insert_inds = get_insert_inds(proteome, peptide_strata[CISSPLICED_KEY])
    for peptide, prot_idx in zip(peptide_strata[CISSPLICED_KEY], insert_inds):
        if prot_idx in modified_ids:
            choices = []
            for i in range(len(proteome[prot_idx])):
                idx_ok = True
                for mod_pos in modified_ids[prot_idx]:
                    if (i-mod_pos) < 20 or (mod_pos-i) < 40:
                        idx_ok = False
                        break
                if idx_ok:
                    choices.append(i)
        else:
            modified_ids[prot_idx] = []
            choices = list(range(len(proteome[prot_idx])))

        proteome[prot_idx], mod_pos, splice_site = add_spliced_seq(
            proteome[prot_idx], peptide, choices
        )
        modified_ids[prot_idx].append(mod_pos)
        peptide_data.append({
            'peptide': peptide,
            'proteinIdx': prot_idx,
            'stratum': CISSPLICED_KEY,
            'frag1': peptide[:splice_site],
            'frag2': peptide[splice_site:],
        })

    insert_inds_1, insert_inds_2 = get_trans_insert_inds(proteome, peptide_strata[TRANSPLICED_KEY])

    for peptide, idx_1, idx_2 in zip(peptide_strata[TRANSPLICED_KEY], insert_inds_1, insert_inds_2):
        mod_idxs = (idx_1, idx_2)
        splice_site = random.randint(2, len(peptide) - 2)
        frags = (peptide[:splice_site], peptide[splice_site:])
        for idx, frag in zip(mod_idxs, frags):
            if idx in modified_ids:
                choices = []
                for i in range(len(proteome[idx])):
                    idx_ok = True
                    for mod_pos in modified_ids[idx]:
                        if i-mod_pos < 20 or mod_pos-i < len(peptide):
                            idx_ok = False
                            break
                    if idx_ok:
                        choices.append(i)
            else:
                modified_ids[idx] = []
                choices = list(range(len(proteome[idx])))
            proteome[idx], mod_pos = add_canonical_seq(
                proteome[idx], frag, choices
            )
            modified_ids[idx].append(mod_pos)
        peptide_data.append({
            'peptide': peptide,
            'proteinIdx': -1,
            'stratum': TRANSPLICED_KEY,
            'frag1': peptide[:splice_site],
            'frag2': peptide[splice_site:],
        })

    peptide_df = pd.DataFrame(peptide_data)

    return proteome, peptide_df, list(modified_ids.keys())
