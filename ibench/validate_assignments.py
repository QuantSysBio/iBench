""" Function for final validation of the assigned strata.
"""
from Bio import SeqIO
import pandas as pd

from ibench.check_presence import (
    find_cis_matched_splice_reactants,
    check_cis_present,
    generate_pairs,
)
from ibench.constants import (
    CANONICAL_KEY,
    CISSPLICED_KEY,
    ENDC_TEXT,
    GT_SCAN_KEY,
    OKCYAN_TEXT,
    TRANSPLICED_KEY,
)

def check_assignment(df_row, modified_proteome, has_cis, enzyme):
    """ Function to check that a peptide is present in the modified proteome according to
        its assigned stratum.

    Parameters
    ----------
    df_row : pd.Series
        Data describing where and how the peptide should be present.
    modified_proteome : list of str
        The proteins in the proteome.
    has_cis : bool
        Flag indicating whether the strata contains cisspliced peptides.
    enzyme : str
        The enzyme used to digest the proteins (either unspecific or trypsin).
    """
    peptide = df_row['peptide']
    if df_row['stratum'] == CANONICAL_KEY:
        protein_idx = df_row['proteinIdx']
        # For discoverable, check that the peptide is present.
        if enzyme == 'trypsin':
            if 'K' + peptide in modified_proteome[protein_idx]:
                return True
        else:
            if peptide in modified_proteome[protein_idx]:
                return True
        return False

    if df_row['stratum'] == CISSPLICED_KEY:
        protein_idx = df_row['proteinIdx']
        # For cisspliced, check that not present as discoverable and is present as cisspliced.
        for protein in modified_proteome:
            if peptide in protein:
                return False
        if check_cis_present(
            modified_proteome[protein_idx],
            df_row['frag1'],
            df_row['frag2'],
        ):
            return True
        return False

    if df_row['stratum'] == TRANSPLICED_KEY:
        # For transspliced, check that not present as discoverable or as cisspliced.
        splice_pairs = generate_pairs(peptide)
        for protein in modified_proteome:
            if peptide in protein:
                return False
            if has_cis and find_cis_matched_splice_reactants(protein, splice_pairs) is not None:
                return False
        return True

    return False


def validate_proteome(hq_df, meta_df, output_folder, enzyme):
    """ Function to run a final validation that all peptides have been properly assigned
        in the proteome and filter any peptides which have not been.

    Parameters
    ----------
    hq_df : pd.DataFrame
        DataFrame of the ground truth peptides.
    meta_df : pd.DataFrame
        DataFrame of meta data related to embedding in the proteome.
    output_fodler : str
        The folder where iBench writes its outputs.
    enzyme : str
        Either unspecific or trypsin.
    """
    print(
        OKCYAN_TEXT +
        f'\tRunning final validation.' +
        ENDC_TEXT
    )
    with open(f'{output_folder}/modified_proteome.fasta', encoding='UTF-8') as prot_file:
        modified_proteome = [
            str(x.seq) for x in SeqIO.parse(prot_file,'fasta')
        ]

    has_cis = CISSPLICED_KEY in meta_df['stratum'].unique().tolist()
    meta_df['valid'] = meta_df.apply(
        lambda x : check_assignment(x, modified_proteome, has_cis, enzyme), axis=1
    )
    n_not_embedded = meta_df[~meta_df['valid']].shape[0]
    print(
        OKCYAN_TEXT +
        f'\tFailed to embed {n_not_embedded} peptides in the proteome.' +
        ENDC_TEXT
    )

    meta_df = meta_df.rename(columns={'peptide': 'il_peptide'})
    hq_df = pd.merge(
        hq_df,
        meta_df[['il_peptide', 'valid']],
        how='inner',
        on='il_peptide',
    )
    hq_df = hq_df[hq_df['valid']]

    hq_df = hq_df.drop('valid', axis=1)
    hq_df = hq_df.sort_values(by=GT_SCAN_KEY)

    hq_df.to_csv(
        f'{output_folder}/high_confidence.csv',
        index=False,
    )

    return n_not_embedded
