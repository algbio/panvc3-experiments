# Copyright (c) Tuukka Norri 2023
# Licenced under the MIT licence.

# vim: syntax=snakefile

from snakemake.utils import min_version
min_version("7.32.4")


# Our cluster does not allow outgoing connections from the nodes, so here are the rules needed to install the environments.
rule bcftools:
	conda:		"../environments/bcftools.yaml"
	shell:		"echo 'Installed dependencies in bcftools.yaml'"


rule biopython:
	conda:		"../environments/biopython.yaml"
	shell:		"echo 'Installed dependencies in biopython.yaml'"


rule bowtie2:
	conda:		"../environments/bowtie2.yaml"
	shell:		"echo 'Installed dependencies in bowtie2.yaml'"


rule bwa:
	conda:		"../environments/bwa.yaml"
	shell:		"echo 'Installed dependencies in bwa.yaml'"


rule gatk:
	conda:		"../environments/gatk.yaml"
	shell:		"echo 'Installed dependencies in gatk.yaml'"


rule graphtyper2:
	conda:		"../environments/graphtyper2.yaml"
	shell:		"echo 'Installed dependencies in graphtyper2.yaml'"


rule gridss:
	conda:		"../environments/gridss.yaml"
	shell:		"echo 'installed dependencies in gridss.yaml'"


rule happy:
	conda:		"../environments/happy.yaml"
	shell:		"echo 'installed dependencies in happy.yaml'"


rule manta:
	conda:		"../environments/manta.yaml"
	shell:		"echo 'installed dependencies in manta.yaml'"


rule mason:
	conda:		"../environments/mason.yaml"
	shell:		"echo 'installed dependencies in mason.yaml'"


rule panvc3:
	conda:		"../environments/panvc3.yaml"
	shell:		"echo 'Installed dependencies in panvc3.yaml'"


rule samtools:
	conda:		"../environments/samtools.yaml"
	shell:		"echo 'Installed dependencies in samtools.yaml'"


rule truvari:
	conda:		"../environments/truvari.yaml"
	shell:		"echo 'Installed dependencies in truvari.yaml'"


rule vg:
	conda:		"../environments/vg.yaml"
	shell:		"echo 'Installed dependencies in vg.yaml'"
