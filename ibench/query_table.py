""" Functions.
"""
from Bio import SeqIO
import pandas as pd
from ibench.check_presence import find_cis_matched_splice_reactants, generate_pairs
from ibench.constants import (
    CANONICAL_KEY,
    CISSPLICED_KEY,
    ENGINE_SCORE_KEY,
    GT_SCAN_KEY,
    GT_SOURCE_KEY,
    LABEL_KEY,
    PEPTIDE_KEY,
    Q_VALUE_KEY,
)
from ibench.input.mascot import read_single_mascot_data
from ibench.input.maxquant import read_single_mq_data
from ibench.input.peaks import read_single_peaks_data
from ibench.input.percolator import read_single_percolator_data

def _remap_to_proteome(peptide, proteome):
    """ Function to check for the presence of an identified peptide as either canonical
        or spliced in the input proteome.
    """
    accession_stratum = 'unknown'
    splice_pairs = generate_pairs(peptide)
    for protein in proteome:
        if peptide in protein:
            return CANONICAL_KEY
        if find_cis_matched_splice_reactants(protein, splice_pairs) is not None:
            accession_stratum = CISSPLICED_KEY
    return accession_stratum

def read_data(location, name, engine, config, flag=''):
    """ Function to read results.
    """
    if engine=='percolator':
        target_df = read_single_percolator_data(
            location,
            1.01,
            -10000,
            hq_hits_only=False,
        )
    elif engine == 'mascot':
        target_df = read_single_mascot_data(
            location, 1.01, -1_000, hq_hits_only=False, filter_ptms=config.filter_ptms,
        )
    elif engine == 'maxquant':
        target_df = read_single_mq_data(
            location, -1_000, hq_hits_only=False, filter_ptms=config.filter_ptms,
        )
    else:
        target_df = read_single_peaks_data(
            location, -1_000, hq_hits_only=False,  filter_ptms=config.filter_ptms,
        )

    if engine != 'percolator':
        if name == 'decoy':
            target_df = target_df[target_df[LABEL_KEY] == -1]
        else:
            target_df = target_df[target_df[LABEL_KEY] == 1]

    target_df = target_df.sort_values(by=ENGINE_SCORE_KEY, ascending=False)
    target_df = target_df.drop_duplicates(subset=['source', 'scan'])
    target_df = target_df.rename(
        columns={
            PEPTIDE_KEY: f'{name}{flag}Peptide',
            ENGINE_SCORE_KEY: f'{name}{flag}Score',
            Q_VALUE_KEY: f'{name}{flag}qValue',
            'source': GT_SOURCE_KEY,
            'scan': GT_SCAN_KEY,
        }
    )
    if config.cisspliced_fraction > 0:
        with open(config.proteome_loc, encoding='UTF-8') as prot_file:
            modified_proteome = [
                str(x.seq) for x in SeqIO.parse(prot_file, 'fasta')
            ]
        target_df[f'{name}{flag}Stratum'] = target_df[f'{name}{flag}Peptide'].apply(
            lambda x : _remap_to_proteome(x, modified_proteome)
        )
    else:
        target_df[f'{name}{flag}Stratum'] = target_df[f'{name}{flag}Peptide'].apply(
            lambda x : CANONICAL_KEY if isinstance(x, str) else None
        )

    output_columns = [
            GT_SOURCE_KEY,
            GT_SCAN_KEY,
            f'{name}{flag}Peptide',
            f'{name}{flag}Score',
            f'{name}{flag}Stratum',
        ]
    if f'{name}qValue' in target_df.columns:
        output_columns.append(f'{name}{flag}qValue')

    target_df = target_df[output_columns]

    return target_df

def create_query_table(config):
    """ Function to create a query table combining ground truth peptides with the identifications
        of one or more methods.

    Parameters
    ----------
    config : ibench.config.Config
        The Config object used to control the experiment.
    """
    qt_df = pd.read_csv(
        f'{config.output_folder}/high_confidence.csv'
    )
    qt_df = qt_df.rename(columns={
        'peptide': 'truePeptide',
        'stratum': 'trueStratum',
    })
    qt_df[GT_SOURCE_KEY] = f'ibenchGroundTruth_{config.identifier}'

    for method in config.benchmark_results:
        name = method['name']
        target_df = read_data(method['resultsLocation'], name, method['searchEngine'], config)

        qt_df = pd.merge(
            qt_df,
            target_df,
            how='left',
            on=[GT_SOURCE_KEY, GT_SCAN_KEY],
        )

        if 'decoyLocation' in method:
            decoy_df = read_data(
                method['decoyLocation'],
                name,
                method['searchEngine'],
                config,
                flag='Decoy',
            )

            qt_df = pd.merge(
                qt_df,
                decoy_df,
                how='left',
                on=[GT_SOURCE_KEY, GT_SCAN_KEY],
            )

    qt_df.to_csv(
        f'{config.output_folder}/queryTable.csv',
        index=False,
    )
