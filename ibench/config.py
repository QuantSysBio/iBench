""" Definition of the Config class which is used to manage the iBench experiement.
"""
import os

import yaml

ALL_CONFIG_KEYS = (
    'benchmarkResults',
    'closenessCutOff',
    'filterPTMs',
    'identifier',
    'scanFolder',
    'scanFormat',
    'searchResults',
    'maxSequenceLength',
    'minSequenceLength',
    'outputFolder',
    'randomSeed',
    'scoreCutOffs',
    'qValueCutOffs',
    'canonicalFraction',
    'cissplicedFraction',
    'transsplicedFraction',
    'proteome',
    'ms2Accuracy',
    'enzyme',
)

class Config:
    """ Holder for configuration of the ibench pipeline.
    """
    def __init__(self, config_file):
        with open(config_file, 'r', encoding='UTF-8') as stream:
            config_dict = yaml.safe_load(stream)
        for config_key in config_dict:
            if config_key not in ALL_CONFIG_KEYS:
                raise ValueError(f'Unrecognised key {config_key} found in config file.')
        self._load_data(config_dict)

    def __str__(self):
        out_str = (
            f'iBench Config for Experiment : {self.identifier}' +
            f'Percentage Canonical : {round(self.discoverable_fraction*100)}' +
            f'Percentage Cisspliced : {round(self.cisspliced_fraction*100)}' +
            f'Percentage Transspliced : {round(self.transspliced_fraction*100)}'
        )
        return out_str

    def _load_data(self, config_dict):
        """ Helper function to load the yaml data provided into the
            Config object.

        Parameters
        ----------
        config_data : dict
            The data parsed from the yaml file.
        """
        self.identifier = config_dict['identifier']
        self.search_results = config_dict.get('searchResults')
        self.output_folder = config_dict.get('outputFolder')
        self.random_seed = config_dict.get('randomSeed', 42)
        self.q_cuts = config_dict.get('qValueCutOffs')
        self.discoverable_fraction = config_dict['canonicalFraction']
        self.cisspliced_fraction = config_dict.get('cissplicedFraction', 0.0)
        self.filter_ptms = config_dict.get('filterPTMs', True)
        self.transspliced_fraction = config_dict['transsplicedFraction']
        self.input_database = config_dict.get('inputDatabase')
        self.proteome_loc = config_dict['proteome']
        self.min_seq_len = config_dict.get('minSequenceLength', 7)
        self.max_seq_len = config_dict.get('maxSequenceLength', 30)
        self.scan_folder = config_dict.get('scanFolder')
        self.scan_format = config_dict.get('scanFormat')
        self.benchmark_results = config_dict.get('benchmarkResults')
        self.ms2_accuracy = config_dict.get('ms2Accuracy', 0.02)
        self.enzyme = config_dict.get('enzyme')

        if self.cisspliced_fraction > 0.0:
            self.closeness_cut_off = config_dict.get('closenessCutOff', 3)
        else:
            self.closeness_cut_off = config_dict.get('closenessCutOff', 1)

        if not os.path.exists(self.output_folder):
            os.makedirs(self.output_folder)

    def validate(self, pipeline):
        """ Check the values in the config file to ensure they are valid.
        """
        if pipeline == 'createDB':
            if not os.path.exists(self.proteome_loc):
                raise ValueError(f'No file at:\n\t{self.proteome_loc}')
        elif pipeline == 'analysis':
            if self.benchmark_results is None:
                raise ValueError('No results provided to benchmark.')
        elif pipeline == 'downloadExample':
            return
        else:
            raise ValueError(f'No such pipeline: {pipeline}')
