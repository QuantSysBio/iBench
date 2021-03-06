""" Function for creating the artificial reference database.
"""
import random
import time

from Bio import SeqIO
import pandas as pd

from ibench.add_seqs import (
    add_sequences,
    add_spliced_seq,
    remove_substring,
)
from ibench.check_presence import (
    check_cis,
    check_cis_present,
    generate_pairs,
)
from ibench.constants import AMINO_ACIDS, CANONICAL_KEY, CISSPLICED_KEY, TRAPPING_KEY
from ibench.utils import get_pepitde_strata
from ibench.validate_assignments import validate_proteome

FASTA_WIDTH = 80

def remove_matches(proteome, peptide_strata, ids_of_interest=None):
    """ Function to remove all canonical matches for all strata from the proteome and
        remove all cisspliced matches for non-canonical strata.

    Parameters
    ----------
    proteome : list of str
        A list of unmodified strata.
    peptide_strata : dict
        A dictionary detailing all canonical, cisspliced, and trapping peptides.
    ids_of_interest : list of int or None
        A list of all the indices to check, if None all indices are considered.

    Returns
    -------
    proteome : list of str
        The input proteome, modified to remove matches.
    modified_ids : list of int
        A list of the indices which were modified.
    """
    modified_ids = []

    if ids_of_interest is None:
        ids_of_interest = set(range(len(proteome)))

    for stratum in peptide_strata:
        for peptide in peptide_strata[stratum]:
            if len(peptide_strata[CISSPLICED_KEY]) and stratum != CANONICAL_KEY:
                splice_pairs = generate_pairs(peptide)

            for idx in ids_of_interest:
                if peptide in proteome[idx]:
                    modified_ids.append(idx)
                    proteome[idx] = remove_substring(proteome[idx], peptide)

                if len(peptide_strata[CISSPLICED_KEY]) and stratum != CANONICAL_KEY:
                    replace_strings = check_cis(proteome[idx], splice_pairs)
                    if replace_strings is not None:
                        for replace_string in replace_strings:
                            modified_ids.append(idx)
                            proteome[idx] = remove_substring(proteome[idx], replace_string)

    return proteome, modified_ids

def modify_db(config):
    """ Function to create the artificial reference database.

    Parameters
    ----------
    config : ibench.config.Config
        The Config object used to manage the experiment.
    """
    hq_df = pd.read_csv(f'{config.output_folder}/high_confidence.csv')
    peptide_strata = get_pepitde_strata(hq_df)

    with open(config.proteome_loc, encoding='UTF-8') as prot_file:
        fasta_sequences = [
            str(x.seq) for x in SeqIO.parse(
                prot_file,
                'fasta'
            )
        ]
    fasta_sequences = [x.replace('I', 'L') for x in fasta_sequences]

    clean_proteome(
        fasta_sequences,
        peptide_strata,
        config.output_folder,
    )

    meta_df = add_seqs_to_proteome(config.output_folder, peptide_strata, config.enzyme)

    validate_proteome(hq_df, meta_df, config.output_folder, config.enzyme)


def check_sequences(
        peptide_df,
        modified_proteome,
        modified_ids,
        enzyme,
        run_idx,
    ):
    """ Function to ensure that all of the sequences are present in the proteome according
        to their assigned strata.
    """
    newly_modified_ids = []

    # Canonical Peptides, check existence and add if not found.
    non_spliced_df = peptide_df[peptide_df['stratum'] == CANONICAL_KEY]
    for _, df_row in non_spliced_df.iterrows():
        protein_idx = df_row['proteinIdx']
        if enzyme == 'trypsin':
            peptide = 'K' + df_row['peptide']
        else:
            peptide = df_row['peptide']
        if run_idx > 1 and protein_idx not in modified_ids:
            continue
        if peptide not in modified_proteome[protein_idx]:
            modified_proteome[protein_idx] += peptide
            newly_modified_ids.append(protein_idx)
            if run_idx > 5:
                print('disc', df_row['peptide'])

    # Cisspliced Peptides, check absence as canonical and presence as spliced.
    cis_spliced_df = peptide_df[peptide_df['stratum'] == CISSPLICED_KEY]
    for _, df_row in cis_spliced_df.iterrows():
        # Check absent as discoverable
        protein_idx = df_row['proteinIdx']
        for other_idx in modified_ids:
            if df_row['peptide'] in modified_proteome[other_idx]:
                if run_idx > 5:
                    print('cis found can', df_row['peptide'])
                modified_proteome[other_idx] = modified_proteome[other_idx].replace(
                    df_row['peptide'],
                    ''.join(
                        [random.choice(AMINO_ACIDS) for _ in range(len(df_row['peptide']))]
                    )
                )
                newly_modified_ids.append(other_idx)

        # If not already validated check existence as a cis-spliced
        if run_idx > 1 and protein_idx not in modified_ids:
            continue

        if not check_cis_present(
                modified_proteome[protein_idx],
                df_row['frag1'],
                df_row['frag2'],
            ):
            if len(df_row['frag1']) > len(df_row['frag2']):
                modified_proteome[protein_idx], _, new_splice_site = add_spliced_seq(
                    modified_proteome[protein_idx],
                    df_row['peptide'],
                    [],
                    splice_site_range=range(1, len(df_row['frag1'])-1)
                )
            else:
                modified_proteome[protein_idx], _, new_splice_site = add_spliced_seq(
                    modified_proteome[protein_idx],
                    df_row['peptide'],
                    [],
                    splice_site_range=range(len(df_row['frag1'])+1, len(df_row['peptide']))
                )
            index = peptide_df.index[peptide_df['peptide'] == df_row['peptide']].tolist()[0]
            peptide_df.loc[index, 'frag1'] = df_row['peptide'][:new_splice_site]
            peptide_df.loc[index, 'frag2'] = df_row['peptide'][new_splice_site:]
            newly_modified_ids.append(protein_idx)
            if run_idx > 5:
                print('cis not found', df_row['peptide'])

    trapping_df = peptide_df[peptide_df['stratum'] == TRAPPING_KEY]
    for _, df_row in trapping_df.iterrows():
        for other_idx in modified_ids:
            # Check existence as discoverable
            if df_row['peptide'] in modified_proteome[other_idx]:
                modified_proteome[other_idx] = modified_proteome[other_idx].replace(
                    df_row['peptide'],
                    ''.join(
                        [random.choice(AMINO_ACIDS) for _ in range(len(df_row['peptide']))]
                    )
                )
                newly_modified_ids.append(other_idx)
                print('Trapping Peptide Present as Canonical:', df_row['peptide'])

            # Check existence as cisspliced
            if cis_spliced_df.shape[0]:
                trans_pairs = generate_pairs(df_row['peptide'])
                replace_frags = check_cis(modified_proteome[other_idx], trans_pairs)
                if replace_frags is not None:
                    for replace_string in replace_frags:
                        newly_modified_ids.append(other_idx)
                        modified_proteome[other_idx] = remove_substring(
                            modified_proteome[other_idx], replace_string
                        )
                    if run_idx > 5:
                        print('Trapping Peptide Present as Spliced:', df_row['peptide'])

    return peptide_df, modified_proteome, newly_modified_ids

def add_seqs_to_proteome(output_folder, peptide_strata, enzyme):
    """ Function to add the ground truth sequences to the proteome.

    Parameters
    ----------
    output_folder : str
        The folder where iBench writes its outputs to.
    peptide_strata : dict
        A dictionary containing all of the ground truth peptides in their strata.
    enzyme : str
        The enzyme used (either nonspecific or trypsin).

    Returns
    -------
    peptide_df : pd.DataFrame
        DataFrame containing meta data on the peptides added to the proteome.
    """
    with open(f'{output_folder}/cleaned_proteome.fasta', encoding='UTF-8') as prot_file:
        cleaned_proteome = [
            str(x.seq) for x in SeqIO.parse(prot_file, 'fasta')
        ]

    modified_proteome, peptide_df, modified_ids = add_sequences(
        cleaned_proteome,
        peptide_strata,
        enzyme,
    )

    idx = 1

    while modified_ids:
        modified_ids = set(modified_ids)
        print(
            f'Sequence Adding, Iteration {idx}, {len(modified_ids)} sequences to check.'
        )
        if idx > 2:
            print(modified_ids)
        peptide_df, modified_proteome, modified_ids = check_sequences(
            peptide_df,
            modified_proteome,
            modified_ids,
            enzyme,
            idx,
        )

        if idx == 30:
            break
        idx += 1

    with open(f'{output_folder}/modified_proteome.fasta', 'w', encoding='UTF-8') as out_file:
        for idx, prot in enumerate(modified_proteome):
            out_file.write(f'>modified_protein_{idx}\n')
            prot_chunks = [prot[i:i+FASTA_WIDTH] for i in range(0, len(prot), FASTA_WIDTH)]
            for chunk in prot_chunks:
                out_file.write(f'{chunk}\n')

    return peptide_df

def clean_proteome(fasta_sequences, peptide_strata, output_folder):
    """ Function for removing all existing matches from the proteome.

    Parameters
    ----------
    fasta_sequences : list of str
        A list of the proteins from which we will remove matches.
    peptide_strata : dict
        A dictionary of the strata and their ground truth peptides.
    output_folder : str
        The folder to which iBench will write its outputs.

    Returns
    -------
    cleaned_proteome : list of str
        The input proteome with all matches removed.
    """
    start = time.time()

    # Clean proteome.
    cleaned_proteome, modified_ids = remove_matches(
        fasta_sequences,
        peptide_strata
    )

    idx = 1
    while modified_ids:
        modified_ids = set(modified_ids)
        print(
            f'Peptide Cleaning, Iteration {idx}, {len(modified_ids)} proteins to clean.'
        )

        cleaned_proteome, modified_ids = remove_matches(
            cleaned_proteome,
            peptide_strata,
            modified_ids,
        )
        idx += 1

    end = time.time()
    print(f'Time spent cleaning proteome: {round(end-start, 2)}s')


    with open(f'{output_folder}/cleaned_proteome.fasta', 'w', encoding='UTF-8') as out_file:
        for idx, prot in enumerate(cleaned_proteome):
            out_file.write(f'>cleaned_protein_{idx}\n')
            out_file.write(f'{prot}\n')

    return cleaned_proteome
