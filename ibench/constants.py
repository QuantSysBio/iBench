""" Constants defined for use throughout the iBench library.
"""
HEADER_TEXT = '\033[95m'
OKBLUE_TEXT = '\033[94m'
OKCYAN_TEXT = '\033[96m'
OKGREEN_TEXT = '\033[92m'
WARNING_TEXT = '\033[93m'
FAIL_TEXT = '\033[91m'
ENDC_TEXT = '\033[0m'
BOLD_TEXT = '\033[1m'
UNDERLINE_TEXT = '\033[4m'

FIGSHARE_PATH = 'https://figshare.com/ndownloader/files/37787349'

ACCESSION_KEY = 'proteins'
CHARGE_KEY = 'charge'
HYDRO_INDEX_KEY = 'hydrophobicityIndex'
MASS_KEY = 'mass'
PEPTIDE_KEY = 'peptide'
INTENSITIES_KEY = 'intensities'
MZS_KEY = 'mzs'
PTM_SEQ_KEY = 'ptm_seq'
SOURCE_KEY = 'source'
SCAN_KEY = 'scan'
LABEL_KEY = 'label'
ENGINE_SCORE_KEY = 'engineScore'
SPECTRAL_ANGLE_KEY = 'spectralAngle'
MASS_DIFF_KEY = 'massDiff'
Q_VALUE_KEY = 'q-value'
DELTA_SCORE_KEY = 'deltaScore'
SEQ_LEN_KEY = 'sequenceLength'
ACCESSION_STRATUM_KEY = 'accessionGroup'
SOURCE_INDEX_KEY = 'sourceIndex'

GT_SCAN_KEY = 'gtScan'
GT_SOURCE_KEY = 'gtSource'
IL_PEPTIDE_KEY = 'il_peptide'

CANONICAL_KEY = 'canonical'
CISSPLICED_KEY = 'cisspliced'
TRANSPLICED_KEY = 'transspliced'

PROTON = 1.007276466622
ELECTRON = 0.00054858
H = 1.007825035
C = 12.0
O = 15.99491463
N = 14.003074

N_TERMINUS = H
C_TERMINUS = O + H
CO = C + O
CHO = C + H + O
NH2 = N + (H * 2)
H2O = (H * 2) + O
NH3 = N + (H * 3)
AMINO_ACIDS = 'ACDEFGHKLMNPQRSTVYW'

ION_OFFSET = {
    'a': N_TERMINUS - CHO,
    'b': N_TERMINUS - H,
    'c': N_TERMINUS + NH2,
    'x': C_TERMINUS + CO - H,
    'y': C_TERMINUS + H,
    'z': C_TERMINUS - NH2,
}

RESIDUE_WEIGHTS = {
    'A': 71.037114,
    'R': 156.101111,
    'N': 114.042927,
    'D': 115.026943,
    'C': 103.009185,
    'E': 129.042593,
    'Q': 128.058578,
    'G': 57.021464,
    'H': 137.058912,
    'I': 113.084064,
    'L': 113.084064,
    'K': 128.094963,
    'M': 131.040485,
    'F': 147.068414,
    'P': 97.052764,
    'S': 87.032028,
    'T': 101.047679,
    'W': 186.079313,
    'Y': 163.06332,
    'V': 99.068414,
}

CONFOUNDING_FEATURE_NAMES = [
    'Peptide Length',
    'Charge',
    'Peptide Mass',
    'Hydrophobicity Index',
    'MS2 Coverage',
    'Signal To Noise'
]
