# iBench
Benchmarking Mass Spectrometry Identification Methods.

iBench is a tool to help you understand the performance of an mass spectrometry identification method

To use iBench you will need:
* A source of ground truth identifications which you are confident in.
* MS data for the ground truth identifications (mgf or mzML format)
* A proteome fasta file to modify.
* An identification method or set of method which you wish to benchmark.

In the example used in this data we have
* Ground Truth Identifications : Synthetic Peptides identified at 1% FDR by both MaxQuant and PEAKS.
* MS data : mgf files for the Synthetic Peptide MS measurments.
* Fasta File : This file comes from the expressed proteome of K562 cell line.
* Identification Method(s) : We provide output files for Mascot search and percolator rescoring for 3 different feature sets.

## Set Up

### Before Downloading

You will also require conda to use the iBench software.

### Setting up your environment:


1) To start with create a new conda environment with python version 3.8:

```
conda create --name ibench python=3.8
```

2) Activate this environment

```
conda activate ibench
```

3) You will then need to install the iBench package:

```
python setup.py install
```

4) To check your installation, run the following command (it is normal for this call to hang for a few seconds on first execution)

```
ibench -h
```

Once you have successfully installed iBench you must run it specifying your pipeline and a config file.

## iBench Execution

iBench is executed with two command line arguments --pipeline which specifies the pipeline to be executed.

### iBench Pipeline: createDB

The first option for the pipeline is "createDB" which processes the ground truth identifications to produce the ground truth datasets, artificial reference database and reindexed MS files needed to benchmark performance against. 

### iBench Pipeline: analysis

After running the first pipeline the user should apply their identification method to the reindexed MS files and artificial reference database. These identifications can be provided to iBench for the "analysis" pipeline. This processes the data and produces a report containing multiple figures describing each method's performance.

### iBench Config File

The second argument provided to iBench is --config_file. This specifies the location of a yaml file which contains all the meta data needed to run the iBench software. Such data includes locations of search results, original MS data, as well as user preferences on the creation of the database and results output. Full details of the configuration settings possible are given in the table at the bottom of this README and an example config file is provided used in the "Running a Small Example" section.

## Running a Small Example

As an example we will demonstrate how iBench could be used to benchmark Mascot with Percolator Rescoring with 3 different feature sets. This example will use ground truth identifications from PEAKS and MaxQuant searches of our synthetic peptide library.


### Create the Ground Truth Dataset

To create the ground truth datasets:

```
ibench --config_file example/config.yml --pipeline createDB
```

### Next Steps

To provide a working example we have performed a Mascot search on the output data and provided

### Running the Analysis

You will now

```
ibench --config_file example/config.yml --pipeline analysis
```

This will create a html report with all relevant plots to compare the performance of the 3 pipelines.

## Full List of iBench Config Settings

This provides all of the configuration setting both optional and required for each pipeline.

### Required for any Execution

iBench will always require an identifier for your experiment and an output folder.

| Key         | Description   |
|---------------|--------------------------|
| identifier    | Some identifier or experiment title for your iBench analysis.   |
| outputFolder  | The folder into which all iBench outputs will be written.   |

### Required for createDB Execution

The following config settings can be used

| Key         | Description                                       |
|-------------|---------------------------------------------------|
| scanFolder  | A folder containing all of the MS files.          |
| scanFormat  | The format of the MS files (either mgf or mzML).  |
| searchResults | A list of all of the search results used. The search results should be a list of outputs from PEAKS, MaxQuant, Mascot, or Percolator, see more details below. |
| canonicalFraction | The fraction of the ground truth identifications which should be embedded as canonical sequences in the artificial reference.   |
| cissplicedFraction | The fraction of the ground truth identifications which should be embedded as cisspliced sequences in the artificial reference.   |
| trappingFraction | The fraction of the ground truth identifications which should be trapping sequences in the artificial reference (not discoverable as canonical or cisspliced).   |
| ms2Accuracy  | The m/z accuracy on the measurement of the MS2 spectrum (needed for calculating coverage and signal to noise features). |
| enzyme       | The enzyme used to produce the peptides, options are unspecific or trypsin. |

### Optional for createDB Execution

The following config settings can be used for specific user requirements in the createDB pipeline.

| Key         | Description                                       |
|-------------|---------------------------------------------------|
| scanFolder  | A folder containing all of the MS files.          |
| randomSeed  | The random seed to use to ensure reproducibility of iBench execution (default=42). |
| maxSequenceLength | Define the maximum sequence length of peptides to be used in the ground truth identifications. |
| minSequenceLength | Define the minimum sequence length of peptides to be used in the ground truth identifications. |
