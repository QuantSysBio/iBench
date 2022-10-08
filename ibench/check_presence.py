""" Utility functions for checking for the presence of a peptide in a proteome.
"""
import re

def find_cis_matched_splice_reactants(protein, splice_reactants):
    """ Function to check if any two potential splice reactants of a peptide are present in
        a single protein.

    Parameters
    ----------
    protein : str
        The protein sequence.
    splice_reactants : list of tuple of str
        A list of the possible splice reactants for the peptide.

    Returns
    -------
    replace_strs : list of str or None
        Either a list of peptide strings which should be replaced in the protein or None
        if the peptide is not present as a cisspliced peptide.
    """
    for (sr_1, sr_2) in splice_reactants:
        replace_strs = []
        if sr_1 in protein and sr_2 in protein:
            frag_1_inds = [
                m.start() for m in re.finditer(sr_1, protein)
            ]
            frag_2_inds = [
                m.start() for m in re.finditer(sr_2, protein)
            ]
            for f1_ind in frag_1_inds:
                for f2_ind in frag_2_inds:
                    if abs(f1_ind - f2_ind) <= 25:
                        if len(sr_1) > 1:
                            replace_strs.append(sr_1)
                        if len(sr_2) > 1:
                            replace_strs.append(sr_2)
                        return replace_strs
    return None

def check_cis_present(protein, sr_1, sr_2):
    """ Function to check if two splice reactants are present in a proteome with an intervening
        sequence length of 25 or less.

    Parameters
    ----------
    protein : str
        The protein sequence.
    sr_1 : str
        The first splice reactant.
    sr_2 : str
        The second splice reactant.

    Returns
    -------
    is_present : bool
        Flag indicating if the two splice reactants are present with an intervening sequence
        length of 25.
    """
    if sr_1 in protein and sr_2 in protein:
        frag_1_inds = [
            m.start() for m in re.finditer(sr_1, protein)
        ]
        frag_2_inds = [
            m.start() for m in re.finditer(sr_2, protein)
        ]
        for f1_ind in frag_1_inds:
            for f2_ind in frag_2_inds:
                if abs(f1_ind - f2_ind) <= 25:
                    return True
    return False


def generate_pairs(peptide):
    """ Function to generate all possible splice reactant pairs for a peptide.

    Parameters
    ----------
    peptide : str
        The peptide for which we wish to generate pairs.

    Returns
    -------
    sr_pairs : list of tuple of str
        A list of the possible splice reactant pairs for the input peptide.
    """
    sr_pairs = [
        (peptide[:i], peptide[i:]) for i in range(1, len(peptide))
    ]
    return sr_pairs
