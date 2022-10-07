""" Utility functions used throughout the iBench library.
"""
import numpy as np

from ibench.constants import (
    CANONICAL_KEY,
    CHARGE_KEY,
    CISSPLICED_KEY,
    ION_OFFSET,
    MZS_KEY,
    PROTON,
    RESIDUE_WEIGHTS,
    TRANSPLICED_KEY,
)

def remove_source_suffixes(source):
    """ Helper function to remove raw, mzML or mgf suffixes from source name.

    Parameters
    ----------
    source : str
        The name of a source file.

    Returns
    -------
    source : str
        The updated name with a suffix removed.
    """
    if source.endswith('.mzML'):
        return source[:-5]
    if source.endswith('.raw') or source.endswith('.mgf'):
        return source[:-4]
    return source

def get_pepitde_strata(hq_df):
    """ Function to create a dictionary of peptide strat from the ground truth DataFrame.

    Parameters
    ----------
    hq_df : pd.DataFrame
        A DataFrame detailing the ground truth peptides.

    Returns
    -------
    peptide_strata : dict
        A dictionary mapping the assigned strata to their peptides.
    """
    peptide_strata = {
        CANONICAL_KEY: hq_df[
            hq_df['stratum'] == CANONICAL_KEY
        ]['il_peptide'].unique().tolist(),
        CISSPLICED_KEY: hq_df[
            hq_df['stratum'] == CISSPLICED_KEY
        ]['il_peptide'].unique().tolist(),
        TRANSPLICED_KEY: hq_df[
            hq_df['stratum'] == TRANSPLICED_KEY
        ]['il_peptide'].unique().tolist(),
    }
    return peptide_strata

def get_matches(df_row, group_masses, frag_z, ms2_accuracy):
    """ Function to get the matched intensities and fragment locations for a set of ions
        of a given charge.

    Parameters
    ----------
    df_row : pd.Series
        DataFrame row containing PSM data.
    group_masses : list of float
        A list of the masses of each possible fragment.
    frag_z : int
        The fragment charge being searched for.
    ms2_accuracy : float
        The accuracy of the mz measurement on the MS2 spectrum.

    Returns
    -------
    matched_locs : list of int
        A list of the matched fragmentation positions.
    """
    matched_locs = []
    matched_inds = []
    for idx, base_mass in enumerate(group_masses):
        fragment_mz = (
            base_mass + (frag_z * PROTON)
        )/frag_z
        matched_mz_ind = np.argmin(
            np.abs(df_row[MZS_KEY] - fragment_mz)
        )
        if np.abs(df_row[MZS_KEY][matched_mz_ind] - fragment_mz) < ms2_accuracy:
            matched_locs.append(idx+1)
            matched_inds.append(matched_mz_ind)

    return matched_locs, matched_inds

def calculate_ms2_feats(df_row, ms2_accuracy):
    """ Function to calculate basic features of the ms2 match for a PSM.

    Parameters
    ----------
    df_row : pd.Series
        DataFrame row containing the PSM spectral data.
    ms2_accuracy : float
        The accuracy of the mz measurement on the MS2 spectrum.

    Returns
    -------
    df_row : pd.Series
        The input row with spectral features calculated.
    """
    # Calculate b,y,a ions
    pep_len = len(df_row['peptide'])
    sub_seq_masses = compute_potential_mzs(
        df_row['peptide'], reverse=False
    )
    rev_sub_seq_masses = compute_potential_mzs(
        df_row['peptide'], reverse=True
    )

    # a- and b-ions for each of the charge options
    possible_ions = {}
    possible_ions['a'] = ION_OFFSET['a'] + sub_seq_masses
    possible_ions['b'] = ION_OFFSET['b'] + sub_seq_masses
    possible_ions['y'] = ION_OFFSET['y'] + rev_sub_seq_masses

    max_charge = min((df_row[CHARGE_KEY], 3))

    # match
    matched_inds = []
    matched_locs = {}
    for ion_group, group_mzs in possible_ions.items():
        matched_locs[ion_group] = []
        for frag_z in range(1, max_charge+1):
            ion_locs, ion_inds = get_matches(
                df_row, group_mzs, frag_z, ms2_accuracy
            )
            if ion_group == 'y':
                matched_locs[ion_group].extend(
                    [pep_len-loc_idx for loc_idx in ion_locs]
                )
            else:
                matched_locs[ion_group].extend(ion_locs)
            matched_inds.extend(ion_inds)

    matched_inds = list(set(matched_inds))

    all_cov_locs = set(matched_locs['a'] + matched_locs['b'] + matched_locs['y'])
    df_row['ms2Coverage'] = len(all_cov_locs)/(pep_len - 1)
    df_row['signalToNoise'] = np.sum(
        np.take(df_row['intensities'], matched_inds)
    )/np.sum(
        df_row['intensities']
    )

    return df_row

def compute_potential_mzs(sequence, reverse):
    """ Function to compute the molecular weights of potential fragments
        generated from a peptide (y & b ions, charges 1,2, or 3, and H2O
        or O2 losses).

    Parameters
    ----------
    sequence : str
        The peptide sequence for which we require molecular weights.
    reverse : bool
        Whether we are getting fragment mzs in the forward direction
        (eg for b ions), or backward direction (eg. for y ions).

    Returns
    -------
    mzs : np.array of floats
        An array of all the possible mzs that coule be observed in
        the MS2 spectrum of a sequence.
    """
    if reverse:
        sequence = sequence[::-1]

    sequence_length = len(sequence)
    n_fragments = sequence_length - 1
    mzs = np.empty(n_fragments)

    tracking_mw = 0.0

    for idx in range(n_fragments):
        tracking_mw += RESIDUE_WEIGHTS[sequence[idx]]
        mzs[idx] = tracking_mw

    tracking_mw += RESIDUE_WEIGHTS[sequence[n_fragments]]

    return mzs
