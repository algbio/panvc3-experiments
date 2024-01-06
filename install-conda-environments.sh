#!/bin/bash

set -euxo pipefail

for target in bcftools biopython bowtie2 bwa gatk gridss happy manta mason panvc2 panvc2_gatk panvc2_vcf2multialign panvc3 samtools truvari vcf2multialign vg
do
	snakemake -c 1 --use-conda --conda-prefix conda-env --snakefile workflow/rules/install.smk -- ${target}
done
