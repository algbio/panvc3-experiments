#!/bin/bash

set -euxo pipefail

for target in bcftools biopython bowtie2 bwa gridss happy manta mason panvc3 samtools truvari vg
do
	snakemake -c 1 --use-conda --conda-prefix conda-env --snakefile workflow/rules/install.smk -- ${target}
done
