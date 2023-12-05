# Copyright (c) Tuukka Norri 2023
# Licenced under the MIT licence.

# vim: syntax=snakefile

# Configuration keys
# ––––––––––––––––––
# alignment_id
# reference (FASTA)
# known_variants (bgzipped VCF)
# chromsomes (list)

rule call:
	message:			"Calling variants with GraphTyper2"
	conda:				"../environments/graphtyper2.yaml"
	threads:			workflow.cores
	input:
		alignments		= f"alignments/{config['alignment_id']}.{{wf}}.bam",
		alignment_index	= f"alignments/{config['alignment_id']}.{{wf}}.bam.bai",
		reference		= config['reference'],
		known_variants	= config['known_variants'],
		regions			= "regions/{regions}.bed"
	output:
		done			= f"graphtyper2/{config['alignment_id']}.{{wf}}.{{regions}}.done"
	benchmark:			f"benchmark/graphtyper2_call/{config['alignment_id']}.{{wf}}.{{regions}}"
	shell:				f"graphtyper genotype {{input.reference}} --sam={{input.alignments}} --prior_vcf={{input.known_variants}} --region_file={{input.regions}} --output=graphtyper2/{config['alignment_id']}.{{wildcards.wf}}.graphtyper2.{{wildcards.regions}} && touch {{output.done}}"


rule list_output:
	message:			"Listing VCF files"
	input:				f"graphtyper2/{config['alignment_id']}.{{wf}}.{{regions}}.done"
	output:				f"graphtyper2/{config['alignment_id']}.{{wf}}.{{regions}}.{{chromosome}}.output-list"
	benchmark:			f"benchmark/graphtyper2_list_output/{config['alignment_id']}.{{wf}}.{{regions}}.{{chromosome}}"
	shell:				f"find graphtyper2/{config['alignment_id']}.{{wildcards.wf}}.{{wildcards.regions}}/{{wildcards.chromosome}} -name '*.vcf.gz' | sort > {{output}}"


rule concatenate_output_lists:
	input:				expand(f"graphtyper2/{config['alignment_id']}.{{{{wf}}}}.{{{{regions}}}}.{{chromosome}}.output-list", chromosome = config['chromosomes'])
	output:				f"graphtyper2/{config['alignment_id']}.{{wf}}.{{regions}}.output-list"
	benchmark:			f"benchmark/graphtyper2_concatenate_output_lists/{config['alignment_id']}.{{wf}}.{{regions}}"
	shell:				"cat {input} > {output}"


rule combine:
	message:			"Concatenating called varaints"
	conda:				"../environments/bcftools.yaml"
	input:				f"graphtyper2/{config['alignment_id']}.{{wf}}.{{regions}}.output-list"
	output:				f"graphtyper2/{config['alignment_id']}.{{wf}}.{{regions}}.bcf"
	benchmark:			f"benchmark/graphtyper2_combine/{config['alignment_id']}.{{wf}}.{{regions}}"
	shell:				"bcftools concat --naive --file-list {input} -Ob -o {output}"


rule normalise_to_vcf_gz:
	message:			"Splitting multiallelics" 
	input:				"{variants}.bcf"
	output:				"{variants}.multiallelics-split.vcf.gz"
	benchmark:			"benchmark/graphtyper2_normalise_to_vcf_gz/{variants}"
	shell:				"bcftools norm --multiallelics - -Oz -o {output} graphtyper2/{input}.bcf"


rule filter:
	message:			"Filtering the called variants"
	conda:				"../environments/biopython.yaml"
	threads:			8
	input:				f"graphtyper2/{config['alignment_id']}.{{wf}}.{{regions}}.multiallelics-split.vcf.gz"
	output:				f"variants/{config['alignment_id']}.{{wf}}.graphtyper2.{{regions}}.vcf.gz"
	benchmark:			f"benchmark/graphtyper2_filter/{config['alignment_id']}.{{wf}}.{{regions}}"
	shell:				"vcf_filter.py --local-script ../workflow/scripts/filter_graphtyper2_output.py {input} AAScore --aa-score 0.5 | bgzip -@ {threads} > {output}"
