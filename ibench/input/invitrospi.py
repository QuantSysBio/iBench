""" Functions for reading inVitroSPI outputs.
"""
import os

from Bio import SeqIO
import pandas as pd
import yaml
from ibench.check_presence import find_cis_matched_splice_reactants, generate_pairs
from ibench.constants import RESIDUE_WEIGHTS
from ibench.input.search_results import combine_dataframes

def get_spectral_angles(config):
    """ Function to run inspire to get spectral angles for invitroSPI identifications.
    """
    inspire_config_data = {
        'experimentTitle': 'Spectral Angle Generation',
        'searchResults': 'ignored.csv',
        'searchEngine': 'mascot',
        'outputFolder': config.output_folder,
        'scansFolder': config.scan_folder,
        'scansFormat': config.scan_format,
        'mzAccuracy': config.inspire_settings['mzAccuracy'],
        'mzUnits': config.inspire_settings['mzUnits'],
        'collisionEnergy': config.inspire_settings['collisionEnergy'],
        'spectralAngleDfs': [f'{config.output_folder}/invitroSPI_temp.csv'],
    }

    with open(
        f'{config.output_folder}/inspire_config.yaml', 'w', encoding='UTF-8'
    ) as inspire_conf_file:
        yaml.dump(inspire_config_data, inspire_conf_file)

    os.system(
        f'inspire --config_file {config.output_folder}/inspire_config.yaml --pipeline spectralAngle'
    )

def fetch_stratum(orig_peptide, product_type, splice_types_str, substrate_prot, all_potential_pcp_mws):
    """ Function to try to check strata of assignments.
    """
    peptide = orig_peptide.replace('I', 'L')
    if product_type.startswith('PCP'):
        return 'canonical'
    splices_types = splice_types_str.split(';')

    if peptide in substrate_prot:
        return 'drop'

    pep_mw = round(sum([
        RESIDUE_WEIGHTS[amino_acid] for amino_acid in peptide
    ]), 1)
    if pep_mw in all_potential_pcp_mws:
        return 'drop'

    for splice_type in splices_types:
        if splice_type == 'cis' or splice_type == 'revCis':
            return 'cisspliced'

    trans_pairs = generate_pairs(peptide)
    cis_matched_srs = find_cis_matched_splice_reactants(
        substrate_prot, trans_pairs, len(substrate_prot)
    )
    if cis_matched_srs is not None:
        return 'drop'

    return 'transspliced'

def generate_all_pcp_mw(substrate_seq):
    """ Function to generate cleavage peptides.
    """
    all_mws = []
    for length in range(5, len(substrate_seq)+1):
        for idx in range(len(substrate_seq)-(length-1)):
            all_mws.append(
                round(sum([
                    RESIDUE_WEIGHTS[substrate_seq[idx+product_idx]] for product_idx in range(length)
                ]), 1)
            )
        print(idx)
        print(substrate_seq[idx:idx+length])
    return set(all_mws)

def read_single_invitrospi_data(results_loc, config, combined_prot):
    """ Function to read in a single ProteasomeDB output file from inVitroSPI.
    """
    prot_db_df = pd.read_csv(f'{results_loc}/ProteasomeDB.csv')
    prot_db_df = prot_db_df[prot_db_df['PTM'].apply(lambda x : not isinstance(x, str))]
    prot_db_df = prot_db_df[prot_db_df['pepSeq'].apply(lambda x : 'C' not in x)]
    prot_db_df['filename'] = prot_db_df['runID'].apply(
        lambda x : x.split('-')[-2]
    )

    sample_df = pd.read_csv(f'{results_loc}/sample_list.csv')
    sample_df = sample_df.rename(
        columns={
            'MSfile': 'source',
            'ionScore': 'engineScore',
        }
    )
    sample_df['source'] = sample_df['source'].apply(
        lambda x : x[:-4]
    )

    relevant_files = os.listdir(config.scan_folder)
    if config.scan_format == 'mgf':
        relevant_files = [x[:-4] for x in relevant_files]
    else:
        relevant_files = [x[:-5] for x in relevant_files]
    sample_df = sample_df[
        sample_df['source'].apply(
            lambda x : x in relevant_files
        )
    ]

    all_id_df = pd.merge(
        prot_db_df,
        sample_df[['filename', 'source']],
        how='inner',
        on='filename',
    )
    substrate_seqs = all_id_df['substrateSeq'].unique().tolist()
    all_potential_pcp_mws = {sub_seq: generate_all_pcp_mw(sub_seq) for sub_seq in substrate_seqs}
    all_id_df = all_id_df.rename(columns={
        'scanNum': 'scan',
        'pepSeq': 'peptide',
    })

    all_id_df['modifiedSequence'] = all_id_df['peptide']

    all_id_df['stratum'] = all_id_df[
        ['peptide', 'productType', 'spliceType', 'substrateSeq']
    ].apply(
        lambda df_row : fetch_stratum(
            df_row['peptide'],
            df_row['productType'],
            df_row['spliceType'],
            combined_prot,
            all_potential_pcp_mws[df_row['substrateSeq']],
        ), axis=1
    )
    all_id_df = all_id_df[all_id_df['stratum'] != 'drop']

    all_id_df.to_csv(
        f'{config.output_folder}/invitroSPI_temp.csv', index=False
    )

    get_spectral_angles(config)

    all_id_df = pd.read_csv(
        f'{config.output_folder}/invitroSPI_temp_spectralAngle.csv'
    )

    os.remove(f'{config.output_folder}/invitroSPI_temp.csv')
    os.remove(f'{config.output_folder}/invitroSPI_temp_spectralAngle.csv')

    all_id_df = all_id_df[all_id_df[['spectralAngle', 'stratum']].apply(
        lambda df_row : df_row['spectralAngle'] > config.spectral_angle_cut_offs[
            df_row['stratum']
        ],
        axis=1
    )]

    return all_id_df

def read_invitro_spi_dbs(config):
    """ Function read invitroSPI results.
    """
    all_results = {}
    combined_prot = ''
    for results in config.search_results:
        results_loc = results['resultsLocation']
        polypep_seqs = [
            str(x.seq) for x in SeqIO.parse(
                f'{results_loc}/polypeptides.fasta', 'fasta'
            )
        ]
        for polypep_seq in polypep_seqs:
            combined_prot += polypep_seq.replace('I', 'L')

    for results in config.search_results:
        search_df = read_single_invitrospi_data(
            results['resultsLocation'],
            config,
            combined_prot,
        )
        if results['identificationGroup'] in all_results:
            all_results[results['identificationGroup']].append(search_df)
        else:
            all_results[results['identificationGroup']] = [search_df]

    return combine_dataframes(all_results, False, config)
