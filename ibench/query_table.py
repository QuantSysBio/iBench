""" Functions.
"""
import multiprocessing as mp

from Bio import SeqIO
import pandas as pd
from ibench.check_presence import (
    find_cis_matched_splice_reactants,
    find_trans_matched,
    generate_pairs,
)
from ibench.constants import (
    CANONICAL_KEY,
    CISSPLICED_KEY,
    ENGINE_SCORE_KEY,
    GT_SCAN_KEY,
    GT_SOURCE_KEY,
    LABEL_KEY,
    PEPTIDE_KEY,
    Q_VALUE_KEY,
    SCAN_KEY,
    SOURCE_KEY,
    TRANSPLICED_KEY,
)
from ibench.input.inspire import read_single_inspire_data
from ibench.input.mascot import read_single_mascot_data
from ibench.input.maxquant import read_single_mq_data
from ibench.input.peaks import read_single_peaks_data
from ibench.input.percolator import read_single_percolator_data

def remap_to_canonical_proteome(peptide, proteome):
    """ Function to remap a peptide to the canonical proteome
    """
    if not isinstance(peptide, str):
        return None
    for protein in proteome:
        if peptide in protein:
            return CANONICAL_KEY
    return TRANSPLICED_KEY

def remap_to_proteome(orig_peptide, proteome, max_intervening, allow_trans=False):
    """ Function to check for the presence of an identified peptide as either canonical
        or spliced in the input proteome.
    """
    accession_stratum = 'unknown'
    peptide = orig_peptide.replace('I', 'L')
    splice_pairs = generate_pairs(peptide.replace('I', 'L'))
    for protein in proteome:
        if peptide in protein:
            return CANONICAL_KEY

    for protein in proteome:
        if find_cis_matched_splice_reactants(protein, splice_pairs, max_intervening) is not None:
            return CISSPLICED_KEY
    if allow_trans:
        if len(proteome) == 1:
            if find_trans_matched(proteome[0], splice_pairs, proteome):
                return TRANSPLICED_KEY
        else:
            return TRANSPLICED_KEY
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
    elif engine == 'inspire':
        target_df = read_single_inspire_data(
            location, 1.01,  -1_000, False,
        )
    else:
        target_df = read_single_peaks_data(
            location, -1_000, hq_hits_only=False,  filter_ptms=config.filter_ptms,
        )
    if engine in ('peaks', 'maxquant'):
        if flag == 'Decoy':
            target_df = target_df[target_df[LABEL_KEY] == -1]
        else:
            target_df = target_df[target_df[LABEL_KEY] == 1]

    if not target_df.shape[0]:
        return None

    target_df = target_df.sort_values(by=ENGINE_SCORE_KEY, ascending=False)
    target_df = target_df.drop_duplicates(subset=['source', 'scan'])
    target_df = target_df.rename(
        columns={
            PEPTIDE_KEY: f'{name}{flag}Peptide',
            ENGINE_SCORE_KEY: f'{name}{flag}Score',
            Q_VALUE_KEY: f'{name}{flag}qValue',
        }
    )

    if target_df['source'].iloc[0].startswith('ibench') or target_df['source'].iloc[0].startswith('controllerType'):
        target_df = target_df.rename(
            columns={
                'source': GT_SOURCE_KEY,
                'scan': GT_SCAN_KEY,
            }
        )
        output_columns = [
            GT_SOURCE_KEY,
            GT_SCAN_KEY,
            f'{name}{flag}Peptide',
            f'{name}{flag}Score',
            f'{name}{flag}Stratum',
        ]
    else:
        output_columns = [
            SOURCE_KEY,
            SCAN_KEY,
            f'{name}{flag}Peptide',
            f'{name}{flag}Score',
            f'{name}{flag}Stratum',
        ]

    if config.cisspliced_fraction > 0 or config.discoverable_fraction == -1:
        with open(f'{config.output_folder}/modified_proteome.fasta', encoding='UTF-8') as prot_file:
            modified_proteome = [
                str(x.seq).replace('I', 'L') for x in SeqIO.parse(prot_file, 'fasta')
            ]

        if flag == 'Decoy':
            modified_proteome = [x[::-1] for x in modified_proteome]

        target_df[f'{name}{flag}Stratum'] = target_df.groupby(f'{name}{flag}Peptide')[f'{name}{flag}Peptide'].transform(
            lambda x : remap_to_proteome(x.iloc[0], modified_proteome, config.max_intervening, allow_trans=config.allow_trans)
        )
        if flag == 'Decoy':
            
            target_df[f'{name}{flag}Stratum'] = target_df[f'{name}{flag}Stratum'].apply(
                lambda x : x if x != 'unknown' else CANONICAL_KEY
            )
    else:
        
        if config.allow_trans:
            with open(f'{config.output_folder}/modified_proteome.fasta', encoding='UTF-8') as prot_file:
                modified_proteome = [
                    str(x.seq).replace('I', 'L') for x in SeqIO.parse(prot_file, 'fasta')
                ]
            target_df[f'{name}{flag}Stratum'] = target_df[f'{name}{flag}Peptide'].apply(
                lambda x : remap_to_canonical_proteome(x, modified_proteome)
            )
        else:
            target_df[f'{name}{flag}Stratum'] = target_df[f'{name}{flag}Peptide'].apply(
                lambda x : CANONICAL_KEY if isinstance(x, str) else None
            )

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

    if config.results_files is not None:
        qt_df = qt_df[qt_df['gtSource'].apply(lambda x : x in config.results_files)]

    qt_df = qt_df.rename(columns={
        'peptide': 'truePeptide',
        'stratum': 'trueStratum',
    })

    scan_files = sorted(qt_df[SOURCE_KEY].unique().tolist())
    if config.single_scan_file:
        qt_df[GT_SOURCE_KEY] = f'ibenchGroundTruth_{config.identifier}'
    else:
        qt_df[GT_SOURCE_KEY] = qt_df[SOURCE_KEY].apply(
            lambda x : f'ibenchGroundTruth_{config.identifier}_{scan_files.index(x)}'
        )


    func_args = [
        (method['resultsLocation'], method['name'], method['searchEngine'], config)
         for method in config.benchmark_results
    ]

    with mp.Pool(processes=10) as pool:
        mp_results = pool.starmap(read_data, func_args)


    for method, target_df in zip(config.benchmark_results, mp_results):
        name = method['name']

        if GT_SOURCE_KEY in target_df.columns:
            merge_keys = [GT_SCAN_KEY]
            target_df = target_df.drop(GT_SOURCE_KEY, axis=1)
        else:
            merge_keys = [SOURCE_KEY, SCAN_KEY]

        qt_df = pd.merge(
            qt_df,
            target_df,
            how='left',
            on=merge_keys,
        )

        if 'decoyLocation' in method:
            decoy_df = read_data(
                method['decoyLocation'],
                name,
                method['searchEngine'],
                config,
                flag='Decoy',
            )
            if decoy_df is not None:
                if GT_SOURCE_KEY in decoy_df.columns:
                    merge_keys = [GT_SCAN_KEY]
                    decoy_df = decoy_df.drop(GT_SOURCE_KEY, axis=1)
                else:
                    merge_keys = [SOURCE_KEY, SCAN_KEY]

                qt_df = pd.merge(
                    qt_df,
                    decoy_df,
                    how='left',
                    on=merge_keys,
                )
        elif method['searchEngine'] in ('peaks', 'maxquant'):
            decoy_df = read_data(
                method['resultsLocation'],
                name,
                method['searchEngine'],
                config,
                flag='Decoy',
            )
            if decoy_df is not None:
                qt_df = pd.merge(
                    qt_df,
                    decoy_df,
                    how='left',
                    on=merge_keys,
                )

    qt_df.to_csv(
        f'{config.output_folder}/queryTable.csv',
        index=False,
    )
