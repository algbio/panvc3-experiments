# Copyright (c) Tuukka Norri 2023
# Licenced under the MIT licence.

# vim: syntax=snakefile

# Configuration keys
# ––––––––––––––––––
# alignment_id
# reference (FASTA)
# known_variants (VCF)
# reads_1 (gzipped FASTQ)
# reads_2 (gzipped FASTQ)
# mem_gb

rule index:
	message:			"Indexing known variants for vg map"
	conda:				"../environments/vg.yaml"
	threads:			workflow.cores
	input:
		reference		= config['reference'],
		known_variants	= f"{config['known_variants']}.gz"
	output:
		index			= multiext("index/vg-map/index", ".gcsa", ".gcsa.lcp", ".xg"),
		temp_dir		= temp(directory("temp/vg"))
	benchmark:			"benchmark/vg_index"
	shell:
		f"mkdir -p temp/vg && vg autoindex --threads {{threads}} --workflow map --tmp-dir {{output.temp_dir}} --target-mem {config['mem_gb']}G -r {{input.reference}} -v {{input.known_variants}} -p index/vg-map/index"
		

rule index_giraffe:
	message:			"Indexing known variants for vg giraffe"
	conda:				"../environments/vg.yaml"
	threads:			workflow.cores
	input:
		reference		= config['reference'],
		known_variants	= f"{config['known_variants']}.gz"
	output:				
		index			= "index/vg-giraffe/index.giraffe.gbz",
		dist			= "index/vg-giraffe/index.dist",
		min			= "index/vg-giraffe/index.min",
		temp_dir		= temp(directory("temp/vg"))
	benchmark:	 		"benchmark/vg_index_giraffe"
	shell:
		f"mkdir -p temp/vg && vg autoindex --threads {{threads}} --workflow giraffe --tmp-dir {{output.temp_dir}} --target-mem {config['mem_gb']}G -r {{input.reference}} -v {{input.known_variants}} -p index/vg-giraffe/index"


rule map:
	conda:			"../environments/vg.yaml"
	threads:		workflow.cores
	input:
		index		= multiext("index/vg-map/index", ".gcsa", ".gcsa.lcp", ".xg"),
		reads_1		= config['reads_1'],
		reads_2		= config['reads_2']
	output:			f"alignments/{config['alignment_id']}.vg-map.bam"
	benchmark:		f"benchmark/vg_map.{config['alignment_id']}"
	shell:
		"vg map --threads {threads} --fastq {input.reads_1} --fastq {input.reads_2} --base-name vg/map/index --surject-to bam > {output}"
	

rule map_giraffe:
	conda:			"../environments/vg.yaml"
	threads:		workflow.cores
	input:
		index		= "index/vg-giraffe/index.giraffe.gbz",
		dist		= "index/vg-giraffe/index.dist",
		min			= "index/vg-giraffe/index.min",
		reads_1		= config['reads_1'],
		reads_2		= config['reads_2']
	output:			f"alignments/{config['alignment_id']}.vg-giraffe.bam"
	benchmark:		f"benchmark/vg_map_giraffe.{config['alignment_id']}"
	shell:
		"vg giraffe --threads {threads} --gbz-name {input.index} --minimizer-name {input.min} --dist-name {input.dist} --fastq-in {input.reads_1} --fastq-in {input.reads_2} -o BAM > {output}"
