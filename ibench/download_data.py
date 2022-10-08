""" Function for downloading the models required for inSPIRE execution.
"""
import os
from urllib.request import urlretrieve
import tarfile

from ibench.constants import ENDC_TEXT, FIGSHARE_PATH, OKCYAN_TEXT


def download_data():
    """ Function to download the example dataset from Figshare

    Parameters
    ----------
    force_reload : bool (default=False)
        Flag indicating whether to remove the existing inSPIRE_models folder
        and redownload all models.
    """

    if os.path.isdir('example'):
        print(
            OKCYAN_TEXT + '\tExample data already downloaded.' + ENDC_TEXT
        )
    else:
        print(
            OKCYAN_TEXT + '\tDownloading data...' + ENDC_TEXT
        )
        urlretrieve(FIGSHARE_PATH, filename=f'{os.getcwd()}/example.tar.gz')
        print(
            OKCYAN_TEXT + '\tExtracting Data...' + ENDC_TEXT
        )
        tar = tarfile.open('example.tar.gz', "r:gz")
        tar.extractall()
        tar.close()
        print(
            OKCYAN_TEXT + '\tDataset ready.' + ENDC_TEXT
        )
