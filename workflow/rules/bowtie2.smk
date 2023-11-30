# Copyright (c) Tuukka Norri 2023
# Licenced under the MIT licence.

# vim: syntax=snakefile

# Configuration keys
# ––––––––––––––––––
# alignment_id
# reference (FASTA)
# reads_1 (gzipped FASTQ)
# reads_2 (gzipped FASTQ)


rule index:
	message:	"Indexing the reference for Bowtie 2"
	conda:		"../environments/bowtie2.yaml"
	benchmark:	"benchmark/bowtie2_index"
	threads:	workflow.cores
	input:		config['reference']
	output:		multiext(f"index/bowtie2/index", ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2")
	shell:		f"bowtie2-build --threads {{threads}} {{input}} index/bowtie2/index"


rule align_reads:
	message:			"Aligning reads with Bowtie 2"
	conda:				"../environments/bowtie2.yaml"
	benchmark:			f"benchmark/bowtie2_align.{config['alignment_id']}"
	threads:			workflow.cores
	input:
		index			= multiext(f"index/bowtie2/index", ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"),
		reads_1			= config['reads_1'],
		reads_2			= config['reads_2']
	output:				f"alignments/{config['alignment_id']}.bowtie2.sam.gz"
	shell:				f"bowtie2 --threads {{threads}} -1 {{input.reads_1}} -2 {{input.reads_2}} -x index/bowtie2/index | gzip > {{output}}"

