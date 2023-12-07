# Copyright (c) Tuukka Norri 2023
# Licenced under the MIT licence.

# vim: syntax=snakefile

rule bgzip:
	conda:		"../environments/samtools.yaml"
	threads:	16
	benchmark:	f"benchmark/bgzip/{{file}}.benchmark"
	input:		"{file}"
	output:		"{file}.gz"
	shell:		"bgzip -k -@ {threads} {input}"


rule bgzip_index:
	conda:		"../environments/samtools.yaml"
	benchmark:	f"benchmark/bgzip_index/{{file}}.benchmark"
	input:		"{file}.gz"
	output:		"{file}.gz.gzi"
	shell:		"bgzip -r {input}"


rule convert_sam_gz_to_bam:
	conda:		"../environments/samtools.yaml"
	threads:	16 # workflow.cores
	benchmark:	f"benchmark/convert_sam_gz_to_bam/{{alignments}}.benchmark"
	input:		"{alignments}.sam.gz"
	output:		"{alignments}.bam"
	shell:		"samtools view -@ {threads} -O BAM -o {output} {input}"
	

rule sort_bam:
	conda:		"../environments/samtools.yaml"
	threads:	16 # workflow.cores
	benchmark:	f"benchmark/sort_bam/{{alignments}}.benchmark"
	input:		"{alignments}.bam"
	output:		"{alignments}.sorted.bam"
	shell:		"../workflow/scripts/set-open-file-limit.sh samtools sort -@ {threads} -o {output} {input}"


rule sort_sam:
	conda:		"../environments/samtools.yaml"
	threads:	16 # workflow.cores
	benchmark:	f"benchmark/sort_sam/{{alignments}}.benchmark"
	input:		"{alignments}.sam"
	output:		"{alignments}.sorted.bam"
	shell:		"../workflow/scripts/set-open-file-limit.sh samtools sort -@ {threads} -o {output} {input}"


rule sort_sam_gz:
	message:	"Sorting the alignments"
	conda:		"../environments/samtools.yaml"
	threads:	16
	benchmark:	f"benchmark/sort_sam_gz/{{alignments}}.benchmark"
	input:		f"{{alignments}}.sam.gz"
	output:		f"{{alignments}}.sorted.bam"
	shell:		"../workflow/scripts/set-open-file-limit.sh samtools sort -@ {threads} -o {output} {input}"


rule sort_by_qname_bam:
	message:	"Sorting the alignments by QNAME"
	conda:		"../environments/samtools.yaml"
	threads:	16
	benchmark:	f"benchmark/sort_by_qname_bam/{{alignments}}.benchmark"
	input:		f"{{alignments}}.bam"
	output:		f"{{alignments}}.qname-sorted.bam"
	shell:		"../workflow/scripts/set-open-file-limit.sh samtools sort -n -@ {threads} -o {output} {input}"


rule sort_by_qname_bam_:
	message:	"Sorting the alignments by QNAME"
	conda:		"../environments/samtools.yaml"
	threads:	16
	benchmark:	f"benchmark/sort_by_qname_bam_/{{alignments}}.benchmark"
	input:		f"{{alignments}}.sorted.bam"
	output:		f"{{alignments}}.qname-sorted.bam"
	shell:		"../workflow/scripts/set-open-file-limit.sh samtools sort -n -@ {threads} -o {output} {input}"


rule sort_by_qname_sam_gz:
	message:	"Sorting the alignments by QNAME"
	conda:		"../environments/samtools.yaml"
	threads:	16
	benchmark:	f"benchmark/sort_by_qname_sam_gz/{{alignments}}.benchmark"
	input:		f"{{alignments}}.sam.gz"
	output:		f"{{alignments}}.qname-sorted.bam"
	shell:		"../workflow/scripts/set-open-file-limit.sh samtools sort -n -@ {threads} -o {output} {input}"


rule index_bam:
	conda:		"../environments/samtools.yaml"
	threads:	16 # workflow.cores
	benchmark:	f"benchmark/index_bam/{{alignments}}.benchmark"
	input:		"{alignments}.bam"
	output:		"{alignments}.bam.bai"
	shell:		"samtools index -@ {threads} {input}"


rule index_fasta_fai:
	conda:		"../environments/samtools.yaml"
	benchmark:	f"benchmark/index_fasta_fai/{{reference}}.benchmark"
	input:		"{reference}.fa"
	output:		"{reference}.fa.fai"
	shell:		"samtools faidx {input}"


rule index_vcf_gz_csi:
	conda:		"../environments/bcftools.yaml"
	benchmark:	f"benchmark/index_vcf_gz_csi/{{variants}}.benchmark"
	input:		"{variants}.vcf.gz"
	output:		"{variants}.vcf.gz.csi"
	shell:		"bcftools index {input}"


rule index_vcf_gz_tbi:
	conda:		"../environments/bcftools.yaml"
	benchmark:	f"benchmark/index_vcf_gz_tbi/{{variants}}.benchmark"
	input:		"{variants}.vcf.gz"
	output:		"{variants}.vcf.gz.tbi"
	shell:		"bcftools index -t {input}"


rule normalise_vcf:
	conda:		"../environments/bcftools.yaml"
	benchmark:	f"benchmark/normalise_vcf/{{variants}}.benchmark"
	input:		"{variants}.vcf"
	output:		"{variants}.normalised.vcf.gz"
	shell:		"bcftools norm -m - -O z -o {output} {input}"


rule normalise_vcf_gz:
	conda:		"../environments/bcftools.yaml"
	benchmark:	f"benchmark/normalise_vcf_gz/{{variants}}.benchmark"
	input:		"{variants}.vcf.gz"
	output:		"{variants}.normalised.vcf.gz"
	shell:		"bcftools norm -m - -O z -o {output} {input}"


ruleorder: sort_bam > sort_sam_gz > sort_sam > sort_by_qname_bam > sort_by_qname_sam_gz > sort_by_qname_bam_
