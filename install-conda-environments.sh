#!/bin/bash

set -euxo pipefail

for target in bcftools biopython bowtie2 bwa gatk graphtyper2 gridss happy manta mason samtools truvari vg
do
	snakemake -c 1 --use-conda --conda-prefix conda-env --snakefile workflow/rules/install.smk -- ${target}
done
