# Copyright (c) Tuukka Norri 2023
# Licenced under the MIT licence.

# vim: syntax=snakefile

# Configuration keys
# ––––––––––––––––––
# reference (FASTA)
# chromosomes (list)
# known_variants_prefix (VCF)
# known_variants_suffix (VCF)
# reads_1 (gzipped FASTQ)
# reads_2 (gzipped FASTQ)
# output_prefix


from snakemake.utils import min_version
min_version("7.32.4")


rule sort_sam_gz:
	message:		"Sorting the alignments"
	conda:			"../environments/samtools.yaml"
	threads:		16
	benchmark:		"{config['output_prefix']}/benchmark/panvc3_sort_sam_gz/{alignments}.benchmark"
	input:			"{alignments}.sam.gz"
	output:			"{alignments}.sorted.bam"
	shell:			"../scripts/set-open-file-limit.sh samtools sort -@ {threads} -o {output} {input}"


rule sort_by_qname_sam_gz:
	message:		"Sorting the alignments by QNAME"
	conda:			"../environments/samtools.yaml"
	threads:		16
	benchmark:		"{config['output_prefix']}/benchmark/panvc3_sort_by_qname_sam_gz/{alignments}.benchmark"
	input:			"{alignments}.sam.gz"
	output:			"{alignments}.qname-sorted.bam"
	shell:			"../scripts/set-open-file-limit.sh samtools sort -n -@ {threads} -o {output} {input}"


rule generate_founder_sequences:
	message:				"Generating founder sequences"
	# FIXME: add Conda
	benchmark:				"{config['output_prefix']}/benchmark/panvc3_vcf2multialign.{chromosome}.f{founder_count}.d{minimum_distance}"
	input:
		reference			= config["reference"],
		variants			= "{config['known_variants_prefix']}{chromosome}{config['known_variants_suffix']}"
	output:
		founders_a2m		= "{config['output_prefix']}/founder-sequences/{chromosome}.f{founder_count}.d{minimum_distance}.a2m.gz"
	shell:
		"vcf2multialign"
		" --founder-sequences={founder_count}"
		" --minimum-distance={minimum_distance}"
		" --input-reference={input.reference}"
		" --reference-sequence={chromosome}"
		" --input-variants={input.variants}"
		" --chromosome={chromosome}"
		" --output-sequences-a2m={output.founders_a2m}"
		" --dst-chromosome={chromosome}"
		" --pipe=run-gzip.sh"


rule filter_reference:
	message:				"Extracting remaining contigs from the reference"
	conda:					"../environments/biopython.yaml"
	benchmark:				"{config['output_prefix']}/benchmark/panvc3_filter_reference"
	threads:				8
	input:
		reference			= config["reference"],
	output:
		remaining_contigs	= "{config['output_prefix']}/founder-sequences/remaining-contigs.fa.gz",
		contig_list			= "{config['output_prefix']}/founder-sequences/contig-list.txt",
	params:
		chromosome_args		= lambda: " ".join(map(lambda x: f"-c {x}", config["chromosomes"]))
	shell:
		"bgzip -c -d -@ {threads} {input.reference} | python3 ../scripts/filter_reference.py {params.chromosome_args} | gzip > {output.remaining_contigs}"


rule combine_indexing_input:
	message:				"Combining reference inputs"
	benchmark:				"{config['output_prefix']}/benchmark/panvc3_combine_indexing_input.f{founder_count}.d{minimum_distance}"
	threads:				8
	input:
		founder_sequences	= expand("{output_prefix}/founder-sequences/{chromosome}.f{{founder_count}}.d{{minimum_distance}}.a2m.gz", output_prefix = config['output_prefix'], chromosome = config['chromosomes']),
		remaining_contigs	= "{config['output_prefix']}/founder-sequences/remaining-contigs.fa.gz"
	output:					
		combined_contigs	= "{config['output_prefix']}/founder-sequences/indexing-input.f{founder_count}.d{minimum_distance}.a2m.gz"
	shell:
		# We don't currently need indexable output, hence we can just concatenate the files. (See also filter_reference.)
		"cat {founder_sequences} {remaininig_contigs} > {output.combined_contigs}"


rule build_msa_index:
	message:				"Building the MSA index"
	conda:					"../environments/panvc3.yaml"
	benchmark:				"{config['output_prefix']}/benchmark/panvc3_index_msa.f{founder_count}.d{minimum_distance}"
	input:					"{config['output_prefix']}/founder-sequences/indexing-input.f{founder_count}.d{minimum_distance}.a2m.gz"
	output:
		index				= "{config['output_prefix']}/msa-index/msa-index.f{founder_count}.d{minimum_distance}.dat",
		unaligned_fasta		= "{config['output_prefix']}/msa-index/unaligned.f{founder_count}.d{minimum_distance}.fa"
	shell:
		"panvc3_index_msa"
		" --build-index"
		" --sequence-inputs={input}"
		" --msa-index-output={output.index}"
		" --output-fasta"
		" --pipe-input='gzip -d -c' > {output.unaligned_fasta}"


rule build_bowtie_index:
	message:	"Indexing the reference for Bowtie 2"
	conda:		"../environments/bowtie2.yaml"
	benchmark:	"{config['output_prefix']}/benchmark/panvc3_index_bowtie2.f{founder_count}.d{minimum_distance}"
	threads:	workflow.cores
	input:		"{config['output_prefix']}/msa-index/unaligned.f{founder_count}.d{minimum_distance}.fa"
	output:		multiext("{config['output_prefix']}/index/bowtie2/index.f{founder_count}.d{minimum_distance}", ".1.bt2l", ".2.bt2l", ".3.bt2l", ".4.bt2l", ".rev.1.bt2l", ".rev.2.bt2l")
	shell:		"bowtie2-build --threads {threads} --large-index {input} {config['output_prefix']}/index/bowtie2/index.f{founder_count}.d{minimum_distance}"


rule bowtie_align_reads:
	message:			"Aligning reads with Bowtie 2"
	conda:				"../environments/bowtie2.yaml"
	benchmark:			"{config['output_prefix']}/benchmark/panvc3_align_bowtie2.f{founder_count}.d{minimum_distance}"
	threads:			workflow.cores
	input:
		index			= multiext("{config['output_prefix']}/index/bowtie2/index.f{founder_count}.d{minimum_distance}", ".1.bt2l", ".2.bt2l", ".3.bt2l", ".4.bt2l", ".rev.1.bt2l", ".rev.2.bt2l"),
		reads_1			= config['reads_1'],
		reads_2			= config['reads_2']
	output:				"alignments/bowtie2/alignments.f{founder_count}.d{minimum_distance}.sam.gz"
	params:
		alignment_count	= lambda wildcards: 2 + wildcards.founder_count # founders + reference + 1
	shell:				"bowtie2 --threads {threads} -k {params.alignment_count} -1 {input.reads_1} -2 {input.reads_2} -x {config['output_prefix']}/index/bowtie2/index.f{founder_count}.d{minimum_distance} | gzip > {output}"


rule project_alignments:
	message:	"Projecting the alignments"
	conda:		"../environments/panvc3.yaml"
	benchmark:	"{config['output_prefix']}/benchmark/panvc3_project_alignments.{aligner}.f{founder_count}.d{minimum_distance}"
	threads:	workflow.cores
	input:		
				reference			= config["reference"],
				msa_index			= "{config['output_prefix']}/msa-index/msa-index.f{founder_count}.d{minimum_distance}.dat",
				seq_output_order	= "{config['output_prefix']}/founder-sequences/contig-list.txt", # FIXME: Use .fai for this.
				alignments			= "{config['output_prefix']}/alignments/{aligner}/alignments.f{founder_count}.d{minimum_distance}.sorted.bam"
	output:		
				alignments			= "{config['output_prefix']}/alignments/{aligner}/alignments.f{founder_count}.d{minimum_distance}.projected.sam.gz"
	shell:		"panvc3_project_alignments"
				" --alignments={input.alignments}"
				" --msa-index={input.msa_index}"
				" --reference={input.reference}"
				" --reference-msa-id=REF"
				" --ref-id-separator=/"
				" --reference-order-input={input.seq_output_order}"
				" --record-index-tag=XI"
				" --preserve-tag=XS"
				" --preserve-tag=YS"
				" | gzip > {output.alignments}"


rule recalculate_mapq:
	message:	"Recalculating MAPQ"
	conda:		"../environments/panvc3.yaml"
	benchmark:	"{config['output_prefix']}/benchmark/panvc3_recalculate_mapq.{aligner}.f{founder_count}.d{minimum_distance}"
	input:		"{config['output_prefix']}/alignments/{aligner}/alignments.f{founder_count}.d{minimum_distance}.projected.qname-sorted.bam"
	output:		"{config['output_prefix']}/alignments/{aligner}/alignments.f{founder_count}.d{minimum_distance}.mapq-recalculated.sam.gz"
	shell:		"panvc3_recalculate_mapq"
				" --alignments={input}"
				" | gzip > {output}"
