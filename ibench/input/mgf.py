""" Functions for reading in scans results in mgf format.
"""
import re

import numpy as np
import pandas as pd
from pyteomics import mgf

from ibench.constants import (
    CHARGE_KEY,
    GT_SCAN_KEY,
    INTENSITIES_KEY,
    MZS_KEY,
    SOURCE_KEY,
)
from ibench.utils import calculate_ms2_feats

def _read_mgf_file(mgf_filename, source, scan_id_mappings, config, file_idx):
    """ Function to process an mgf file to find matches with selected scan IDs
        and write the iBench ground truth mgf file.

    Parameters
    ----------
    mgf_filename : str
        The mgf file from which we are reading.
    scan_ids : list of int
        A list of the scan IDs we require.
    scan_file_format : str
        The format of the file used.
    source_list : list of str
        A list of source names.

    Returns
    -------
    scans_df : pd.DataFrame
        A DataFrame of scan results.
    """
    matched_scan_ids = []
    matched_intensities = []
    matched_mzs = []
    precursor_charges = []

    new_spectra = []
    with mgf.read(mgf_filename) as reader:
        for spectrum in reader:
            regex_match = re.match(
                r'(\d+)(.*?)',
                spectrum['params']['title'].split('scan=')[-1]
            )
            scan_id = int(regex_match.group(1))

            if scan_id in scan_id_mappings:
                spectrum['params']['title'] = spectrum['params']['title'].replace(
                    f'scan={scan_id}', f'scan={scan_id_mappings[scan_id]}'
                ).replace(
                    source, f'ibenchGroundTruth_{config.identifier}'
                )
                spectrum['params']['scans'] = scan_id_mappings[scan_id]
                new_spectra.append(spectrum)
                matched_scan_ids.append(scan_id_mappings[scan_id])
                precursor_charges.append(spectrum['params']['charge'][0])
                matched_intensities.append(np.array(list(spectrum['intensity array'])))
                matched_mzs.append(np.array(list(spectrum['m/z array'])))

    if file_idx == 0:
        file_mode = 'w'
    else:
        file_mode = 'a'

    mgf.write(
        new_spectra,
        output=f'{config.output_folder}/ibenchGroundTruth_{config.identifier}.mgf',
        file_mode=file_mode,
    )

    mgf_df = pd.DataFrame(
        {
            GT_SCAN_KEY: pd.Series(matched_scan_ids),
            CHARGE_KEY: pd.Series(precursor_charges),
            INTENSITIES_KEY: pd.Series(matched_intensities),
            MZS_KEY: pd.Series(matched_mzs)
        }
    )

    mgf_df = mgf_df.drop_duplicates(
        subset=[GT_SCAN_KEY]
    )

    return mgf_df

def process_mgf_files(hq_df, config):
    """ Function to read in mgf files, combine ms2 spectral data with the ground truth
        dataset, and write a reindexed mgf file.

    Parameters
    ----------
    hq_df : pd.DataFrame
        The DataFrame of high quality PSMs identified in the original search.
    config : ibench.config.Config
        The Config object which controls the experiment.

    Returns
    -------
    combined_df : pd.DataFrame
        The input DataFrame with additional information gathered from each mgf file.
    """
    sub_df_list = []
    source_files = hq_df[SOURCE_KEY].unique().tolist()
    for file_idx, source_name in enumerate(source_files):
        mgf_file = f'{config.scan_folder}/{source_name}.mgf'
        sub_df = hq_df[hq_df['source'] == source_name]
        scan_mappings = dict(zip(sub_df['scan'].tolist(), sub_df[GT_SCAN_KEY].tolist()))
        mgf_df = _read_mgf_file(mgf_file, source_name, scan_mappings, config, file_idx)

        sub_df = pd.merge(
            sub_df,
            mgf_df,
            how='inner',
            on=GT_SCAN_KEY,
        )

        if sub_df.shape[0]:
            sub_df = sub_df.apply(
                lambda x : calculate_ms2_feats(x, config.ms2_accuracy),
                axis=1,
            )
            sub_df = sub_df.drop([MZS_KEY, 'intensities'], axis=1)
            sub_df_list.append(sub_df)

    return pd.concat(sub_df_list)
