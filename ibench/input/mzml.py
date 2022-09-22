""" Functions for loading experimental spectra from mzml files.
"""
import re

import numpy as np
import pandas as pd
from pyteomics import mzml

from ibench.constants import (
    CHARGE_KEY,
    INTENSITIES_KEY,
    MZS_KEY,
    GT_SCAN_KEY,
    SCAN_KEY,
    SOURCE_KEY,
)
from ibench.utils import calculate_ms2_feats

def _read_mzml_file(mzml_filename, source, scan_id_mappings, config, file_idx):
    """ Function to process an mzml file to find matches with selected scan IDs
        and write the iBench ground truth mzml file.

    Parameters
    ----------
    mzml_filename : str
        The mzml file from which we are reading.
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
    with mzml.read(mzml_filename) as reader:
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
                new_spectra.append(spectrum)
                matched_scan_ids.append(scan_id_mappings[scan_id])
                precursor_charges.append(spectrum['params']['charge'][0])
                matched_intensities.append(np.array(list(spectrum['intensity array'])))
                matched_mzs.append(np.array(list(spectrum['m/z array'])))

    mzml_df = pd.DataFrame(
        {
            GT_SCAN_KEY: pd.Series(matched_scan_ids),
            CHARGE_KEY: pd.Series(precursor_charges),
            INTENSITIES_KEY: pd.Series(matched_intensities),
            MZS_KEY: pd.Series(matched_mzs)
        }
    )

    mzml_df = mzml_df.drop_duplicates(
        subset=[GT_SCAN_KEY]
    )

    return mzml_df, new_spectra

def process_mzml_files(hq_df, config):
    """ Function to read in mzml files, combine ms2 spectral data with the ground truth
        dataset, and write a reindexed mzml file.

    Parameters
    ----------
    hq_df : pd.DataFrame
        The DataFrame of high quality PSMs identified in the original search.
    config : ibench.config.Config
        The Config object which controls the experiment.

    Returns
    -------
    combined_df : pd.DataFrame
        The input DataFrame with additional information gathered from each mzml file.
    """
    sub_df_list = []
    all_spectra = []
    source_files = hq_df[SOURCE_KEY].unique().tolist()
    for file_idx, source_name in enumerate(source_files):
        mzml_file = f'{config.scan_folder}/{source_name}.mzML'
        source_name = mzml_file.split('/')[-1].strip('.mzML')
        sub_df = hq_df[hq_df['source'] == source_name]
        scan_mappings = dict(zip(sub_df['scan'].tolist(), sub_df[GT_SCAN_KEY].tolist()))
        mzml_df, new_spectra = _read_mzml_file(mzml_file, source_name, scan_mappings, config, file_idx)
        all_spectra.extend(new_spectra)
        sub_df = pd.merge(
            sub_df,
            mzml_df,
            how='inner',
            on=GT_SCAN_KEY,
        )

        if sub_df.shape[0]:
            sub_df = sub_df.apply(
                lambda x : calculate_ms2_feats(x, config.ms2_accuracy),
                axis=1,
            )
            sub_df = sub_df.drop(['mzs', 'intensities'], axis=1)
            sub_df_list.append(sub_df)

    mzml.write(
        all_spectra,
        output=f'{config.output_folder}/ibenchGroundTruth_{config.identifier}.mzML',
        file_mode='w',
    )

    return pd.concat(sub_df_list)
