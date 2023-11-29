# Copyright (c) Tuukka Norri 2023
# Licenced under the MIT licence.

# vim: syntax=snakefile

# Configuration keys
# ––––––––––––––––––
# alignment_id
# reference (FASTA)
# chromosomes (list)
# known_variants_prefix (VCF)
# known_variants_suffix (VCF)
# reads_1 (gzipped FASTQ)
# reads_2 (gzipped FASTQ)
# panvc3_conda_environment_path
# vcf2multialign_conda_environment_path

# No global variables, see https://github.com/snakemake/snakemake/issues/2499


from snakemake.utils import min_version
min_version("7.32.4")


def none_if_empty(x):
	return x if bool(x) else None


wildcard_constraints:
	founder_count		= r"\d+",
	minimum_distance	= r"\d+"


rule sort_sam_gz:
	message:		"Sorting the alignments"
	conda:			"../environments/samtools.yaml"
	threads:		16
	benchmark:		f"benchmark/panvc3/sort_sam_gz/{{alignments}}.benchmark"
	input:			f"{{alignments}}.sam.gz"
	output:			f"{{alignments}}.sorted.bam"
	shell:			"../scripts/set-open-file-limit.sh samtools sort -@ {threads} -o {output} {input}"


rule sort_by_qname_sam_gz:
	message:		"Sorting the alignments by QNAME"
	conda:			"../environments/samtools.yaml"
	threads:		16
	benchmark:		f"benchmark/panvc3/sort_by_qname_sam_gz/{{alignments}}.benchmark"
	input:			f"{{alignments}}.sam.gz"
	output:			f"{{alignments}}.qname-sorted.bam"
	shell:			"../scripts/set-open-file-limit.sh samtools sort -n -@ {threads} -o {output} {input}"


rule generate_founder_sequences:
	message:				"Generating founder sequences"
	conda:					none_if_empty(config['vcf2multialign_conda_environment_path'])
	benchmark:				f"benchmark/panvc3_vcf2multialign.{{chromosome}}.f{{founder_count}}.d{{minimum_distance}}"
	input:
		reference			= config["reference"],
		variants			= f"{config['known_variants_prefix']}{{chromosome}}{config['known_variants_suffix']}"
	output:
		founders_a2m		= f"panvc3/founder-sequences/chromosome.{{chromosome}}.f{{founder_count}}.d{{minimum_distance}}.a2m.gz"
	shell:
		"vcf2multialign"
		" --founder-sequences={wildcards.founder_count}"
		" --minimum-distance={wildcards.minimum_distance}"
		" --input-reference={input.reference}"
		" --reference-sequence={wildcards.chromosome}"
		" --input-variants={input.variants}"
		" --chromosome={wildcards.chromosome}"
		" --output-sequences-a2m={output.founders_a2m}"
		" --dst-chromosome={wildcards.chromosome}"
		" --pipe=../workflow/scripts/run-gzip.sh"


rule filter_reference:
	message:				"Extracting remaining contigs from the reference"
	conda:					"../environments/biopython.yaml"
	benchmark:				f"benchmark/panvc3/filter_reference"
	threads:				8
	input:
		reference			= config["reference"],
	output:
		remaining_contigs	= f"panvc3/founder-sequences/remaining-contigs.fa.gz",
		contig_list			= f"panvc3/founder-sequences/contig-list.txt",
	params:
		chromosome_args		= lambda _: " ".join(map(lambda x: f"-c {x}", config["chromosomes"]))
	shell:
		"bgzip -c -d -@ {threads} {input.reference} | python3 ../scripts/filter_reference.py {params.chromosome_args} | gzip > {output.remaining_contigs}"


rule combine_indexing_input:
	message:				"Combining reference inputs"
	benchmark:				f"benchmark/panvc3/combine_indexing_input.f{{founder_count}}.d{{minimum_distance}}"
	threads:				8
	input:
		founder_sequences	= expand("panvc3/founder-sequences/chromosome.{chromosome}.f{{founder_count}}.d{{minimum_distance}}.a2m.gz", chromosome = config['chromosomes']),
		remaining_contigs	= f"panvc3/founder-sequences/remaining-contigs.fa.gz"
	output:					
		combined_contigs	= f"panvc3/founder-sequences/indexing-input.f{{founder_count}}.d{{minimum_distance}}.a2m.gz"
	shell:
		# We don't currently need indexable output, hence we can just concatenate the files. (See also filter_reference.)
		"cat {input.founder_sequences} {input.remaining_contigs} > {output.combined_contigs}"


rule build_msa_index:
	message:				"Building the MSA index"
	conda:					none_if_empty(config['panvc3_conda_environment_path'])
	benchmark:				f"benchmark/panvc3/index_msa.f{{founder_count}}.d{{minimum_distance}}"
	input:					f"panvc3/founder-sequences/indexing-input.f{{founder_count}}.d{{minimum_distance}}.a2m.gz"
	output:
		index				= f"panvc3/msa-index/msa-index.f{{founder_count}}.d{{minimum_distance}}.dat",
		unaligned_fasta		= f"panvc3/msa-index/unaligned.f{{founder_count}}.d{{minimum_distance}}.fa"
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
	benchmark:	f"benchmark/panvc3/index_bowtie2.f{{founder_count}}.d{{minimum_distance}}"
	threads:	workflow.cores
	input:		f"panvc3/msa-index/unaligned.f{{founder_count}}.d{{minimum_distance}}.fa"
	output:		multiext(f"index/panvc3/bowtie2/index.f{{founder_count}}.d{{minimum_distance}}", ".1.bt2l", ".2.bt2l", ".3.bt2l", ".4.bt2l", ".rev.1.bt2l", ".rev.2.bt2l")
	shell:		f"bowtie2-build --threads {{threads}} --large-index {{input}} index/panvc3/bowtie2/index.f{{wildcards.founder_count}}.d{{wildcards.minimum_distance}}"


rule bowtie_align_reads:
	message:			"Aligning reads with Bowtie 2"
	conda:				"../environments/bowtie2.yaml"
	benchmark:			f"benchmark/panvc3/align.{config['alignment_id']}.bowtie2.f{{founder_count}}.d{{minimum_distance}}"
	threads:			workflow.cores
	input:
		index			= multiext(f"index/panvc3/bowtie2/index.f{{founder_count}}.d{{minimum_distance}}", ".1.bt2l", ".2.bt2l", ".3.bt2l", ".4.bt2l", ".rev.1.bt2l", ".rev.2.bt2l"),
		reads_1			= config['reads_1'],
		reads_2			= config['reads_2']
	output:				f"alignments/{config['alignment_id']}.panvc3-bowtie2-f{{founder_count}}-d{{minimum_distance}}.sam.gz"
	params:
		alignment_count	= lambda wildcards: 2 + int(wildcards.founder_count) # founders + reference + 1
	shell:				f"bowtie2 --threads {{threads}} -k {{params.alignment_count}} -1 {{input.reads_1}} -2 {{input.reads_2}} -x index/panvc3/bowtie2/index.f{{wildcards.founder_count}}.d{{wildcards.minimum_distance}} | gzip > {{output}}"


rule project_alignments:
	message:	"Projecting the alignments"
	conda:		none_if_empty(config['panvc3_conda_environment_path'])
	benchmark:	f"benchmark/panvc3/project_alignments.{config['alignment_id']}.{{aligner}}.f{{founder_count}}.d{{minimum_distance}}"
	threads:	workflow.cores
	input:		
				reference			= config["reference"],
				msa_index			= f"panvc3/msa-index/msa-index.f{{founder_count}}.d{{minimum_distance}}.dat",
				seq_output_order	= f"panvc3/founder-sequences/contig-list.txt", # FIXME: Use .fai for this.
				alignments			= f"alignments/{config['alignment_id']}.panvc3-{{aligner}}-f{{founder_count}}-d{{minimum_distance}}.sorted.bam"
	output:		
				alignments			= f"alignments/{config['alignment_id']}.panvc3-{{aligner}}-f{{founder_count}}-d{{minimum_distance}}.projected.sam.gz"
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
	conda:		none_if_empty(config['panvc3_conda_environment_path'])
	threads:	3
	benchmark:	f"benchmark/panvc3/recalculate_mapq.{config['alignment_id']}.{{aligner}}.f{{founder_count}}.d{{minimum_distance}}"
	input:		f"alignments/{config['alignment_id']}.panvc3-{{aligner}}-f{{founder_count}}-d{{minimum_distance}}.projected.qname-sorted.bam"
	output:		f"alignments/{config['alignment_id']}.panvc3-{{aligner}}-f{{founder_count}}-d{{minimum_distance}}.mapq-recalculated.sam.gz"
	shell:		"panvc3_recalculate_mapq"
				" --alignments={input}"
				" | gzip > {output}"


rule max_mapq:
	message:	"Filtering alignments by maximum MAPQ"
	conda:		none_if_empty(config['panvc3_conda_environment_path'])
	threads:	3
	benchmark:	f"benchmark/panvc3/max_mapq.{config['alignment_id']}.{{aligner}}.f{{founder_count}}.d{{minimum_distance}}"
	input:		f"alignments/{config['alignment_id']}.panvc3-{{aligner}}-f{{founder_count}}-d{{minimum_distance}}.mapq-recalculated.sam.gz"
	output:		f"alignments/{config['alignment_id']}.panvc3-{{aligner}}-f{{founder_count}}-d{{minimum_distance}}.max-mapq.sam.gz"
	shell:		"panvc3_subset_alignments"
				" --alignments={input}"
				" --best-mapq"
				" | gzip > {output}"


# Rewrite the CIGAR strings for e.g. Manta.
rule alignment_match:
	message:	"Rewriting CIGAR strings to use alignment match operations"
	conda:		none_if_empty(config['panvc3_conda_environment_path'])
	threads:	3
	benchmark:  f"benchmark/panvc3/alignment_match/{{alignments}}.benchmark"
	input:		"{alignments}.sam.gz"
	output:		"{alignments}.alignment-match.sam.gz"
	shell:		"panvc3_rewrite_cigar"
				" --alignments={input}"
				" --output-alignment-match-ops"
				" | gzip > {output}"
