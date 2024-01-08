# [PanVC 3](https://github.com/tsnorri/panvc3) Experiment Scripts

> **_NOTE:_** This repository contains the experiments described in our reviesed manuscript completed in the beginning of January 2024. In case there are any problems running the experiments, please contact us by e-mail.

## Requirements
 
 * [Snakemake](https://snakemake.readthedocs.io/) (tested with version 7.32.4 and Python 3.11)
 * [Conda](https://conda.io/) (tested with version 22.11.1)

> **_NOTE:_** Using the experiment scripts with Python 3.12 may require Snakemake 8.1 due to [a change in f-string handling](https://github.com/snakemake/snakemake/issues/2586).

## Preparing the experiments

1. Clone the repository with `git clone --recursive https://github.com/algbio/panvc3-experiments.git`
2. Run `./setup.sh` to download the experiment data (including the hs37d5 reference and the phase 3 variant data from The 1000 Genomes Project)
3. Run `cd lib/reference_flow && sh src/download_1kg_pop_table.sh` to download the data files required by Reference Flow.

## Running the experiments

Each of the following directories contains a Snakefile that runs the experiment, as well as a script called `run-snakemake-test.sh`. The script tests the experiment workflow by passing recommended parameters and `--dry-run` to Snakemake. Removing `--dry-run` results in a command suitable for running each experiment.

- [reference-bias](reference-bias)
- [svcalls](svcalls)
- [take-one-out-chr1](take-one-out-chr1)

In all cases, the default target of the Snakefile summarises the results. The Snakefiles contain some variables specific to the testing environment (such as `MEM_GB` set to 168 to specify the memory limit in gigabytes for vg) and hence may have to be modified before running the experiment.

The repository contains Conda environment specifications such that the required tools are installed automatically when running Snakemake.
