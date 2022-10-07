""" Functions for loading experimental spectra from mzml files.
"""
import numpy as np
import pandas as pd
from pyopenms import MSExperiment, MzMLFile # pylint: disable-msg=E0611

from ibench.constants import (
    CHARGE_KEY,
    INTENSITIES_KEY,
    MZS_KEY,
    GT_SCAN_KEY,
    SCAN_KEY,
    SOURCE_KEY,
)
from ibench.utils import calculate_ms2_feats

def _read_mzml_file(mzml_filename, scan_id_mappings, new_exp):
    """ Function to process an MzML file to find matches with scan IDs.
    Parameters
    ----------
    mzml_filename : str
        The mzml file from which we are reading.
    scan_ids : list of int
        A list of the scan IDs we require.
    Returns
    -------
    scans_df : pd.DataFrame
        A DataFrame of scan results.
    """
    matched_scan_ids = []
    matched_intensities = []
    matched_mzs = []
    precursor_charges = []

    exp = MSExperiment()
    MzMLFile().load(mzml_filename, exp)
    for spectrum in exp:
        scan_id = int(spectrum.getNativeID().split('scan=')[1])

        if scan_id in scan_id_mappings:
            replacement_native_id = spectrum.getNativeID().replace(
                f'scan={scan_id}', f'scan={scan_id_mappings[scan_id]}'
            )
            spectrum.setNativeID(replacement_native_id)
            new_exp.addSpectrum(spectrum)
            matched_scan_ids.append(scan_id_mappings[scan_id])
            matched_intensities.append(np.array([peak.getIntensity() for peak in spectrum]))
            matched_mzs.append(np.array([peak.getMZ() for peak in spectrum]))
            precursor_charges.append(spectrum.getPrecursors()[0].getCharge())


    scans_df =  pd.DataFrame(
        {
            GT_SCAN_KEY: pd.Series(matched_scan_ids),
            CHARGE_KEY: pd.Series(precursor_charges),
            INTENSITIES_KEY: pd.Series(matched_intensities),
            MZS_KEY: pd.Series(matched_mzs)
        }
    )

    scans_df = scans_df.drop_duplicates(subset=[GT_SCAN_KEY])

    return scans_df, new_exp


def process_mzml_files(hq_df, config):
    """ Function to read in mzML files, combine ms2 spectral data with the ground truth
        dataset, and write a reindexed mzML file.
    Parameters
    ----------
    hq_df : pd.DataFrame
        The DataFrame of high quality PSMs identified in the original search.
    config : ibench.config.Config
        The Config object which controls the experiment.
    Returns
    -------
    combined_df : pd.DataFrame
        The input DataFrame with additional information gathered from each mzML file.
    """
    sub_df_list = []
    mzml_exp = MSExperiment()
    source_files = hq_df[SOURCE_KEY].unique().tolist()
    for source_name in source_files:
        mzml_file = f'{config.scan_folder}/{source_name}.mzML'

        sub_df = hq_df[hq_df[SOURCE_KEY] == source_name]
        scan_mappings = dict(zip(sub_df[SCAN_KEY].tolist(), sub_df[GT_SCAN_KEY].tolist()))
        mzml_df, mzml_exp = _read_mzml_file(mzml_file, scan_mappings, mzml_exp)

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

    MzMLFile().store(
        f'{config.output_folder}/ibenchGroundTruth_{config.identifier}.mzML',
        mzml_exp
    )

    with open(
        f'{config.output_folder}/ibenchGroundTruth_{config.identifier}.mzML',
        'r',
        encoding='UTF-8',
    ) as file :
        file_data = file.read()

    source_files = hq_df[SOURCE_KEY].unique().tolist()
    for source_name in source_files:
        mzml_file = f'{config.scan_folder}/{source_name}.mzML'
        file_data = file_data.replace(
            source_name, f'ibenchGroundTruth_{config.identifier}'
        )

    with open(
        f'{config.output_folder}/ibenchGroundTruth_{config.identifier}.mzML',
        'w',
        encoding='UTF-8',
    ) as file:
        file.write(file_data)

    return pd.concat(sub_df_list)
