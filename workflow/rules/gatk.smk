# Copyright (c) Tuukka Norri 2023
# Licenced under the MIT licence.

# vim: syntax=snakefile

# Configuration keys
# ––––––––––––––––––
# alignment_id
# reference (FASTA)
# known_variants (bgzipped VCF)
# chromsomes (list)
# ploidy
# reads_1
# reads_2
# mem_mb

rule index_ref:
	message:			"Creating sequence dictionary"
	conda:				"../environments/gatk.yaml"
	input:				"{reference}.fa"
	output:				"{reference}.dict"
	benchmark:			"benchmark/gatk_index_ref/{reference}"
	shell:				"gatk"
						f" --java-options '-Xmx{config['mem_mb']}M -Dsnappy.disable=true -Dsamjdk.snappy.disable=true'"
						" CreateSequenceDictionary"
						" --REFERENCE {input}"
						" --OUTPUT {output}"


rule sort_coordinate:
	message:			"Sorting the alignments"
	conda:				"../environments/gatk.yaml"
	input:				"{alignments}.bam"
	benchmark:			"benchmark/gatk_sort_coordinate/{alignments}"
	output:
		alignments		= "{alignments}.gatk-sorted-coordinate.bam",
		tempdir			= temp(directory("temp/gatk/sort_coordinate/{alignments}"))
	shell:				"mkdir -p {output.tempdir} && gatk"
						f" --java-options '-Xmx{config['mem_mb']}M -Djava.io.tmpdir={{output.tempdir}} -Dsnappy.disable=true -Dsamjdk.snappy.disable=true'"
						" SortSam"
						" --TMP_DIR {output.tempdir}"
						" --INPUT {input}"
						" --OUTPUT {output.alignments}"
						" --SORT_ORDER coordinate"
						" --CREATE_INDEX false"
						" --CREATE_MD5_FILE false"
						" --USE_JDK_DEFLATER true"
						" --USE_JDK_INFLATER true"


rule sort_queryname:
	message:			"Sorting the alignments"
	conda:				"../environments/gatk.yaml"
	input:				"{alignments}.bam"
	output:
		alignments		= "{alignments}.gatk-qname-sorted.bam",
		tempdir			= temp(directory("temp/gatk/sort_queryname/{alignments}"))
	benchmark:			"benchmark/gatk_sort_queryname/{alignments}"
	shell:				"mkdir -p {output.tempdir} && gatk"
						f" --java-options '-Xmx{config['mem_mb']}M -Djava.io.tmpdir={{output.tempdir}} -Dsnappy.disable=true -Dsamjdk.snappy.disable=true'"
						" SortSam"
						" --TMP_DIR {output.tempdir}"
						" --INPUT {input}"
						" --OUTPUT {output.alignments)"
						" --SORT_ORDER queryname"
						" --CREATE_INDEX false"
						" --CREATE_MD5_FILE false"
						" --USE_JDK_DEFLATER true"
						" --USE_JDK_INFLATER true"


rule fastq_to_unaligned_bam:
	message:			"Converting FASTQ to unaligned BAM"
	conda:				"../environments/gatk.yaml"
	input:
		reads_1			= config["reads_1"],
		reads_2			= config["reads_2"]
	output:				
		alignments		= f"gatk/{config['alignment_id']}.unaligned.bam",
		tempdir			= temp(directory(f"temp/gatk/fastq_to_unaligned_bam/{config['alignment_id']}"))
	benchmark:			f"benchmark/gatk_fastq_to_unaligned_bam/{config['alignment_id']}"
	shell:				"mkdir -p {output.tempdir} && gatk"
						f" --java-options '-Xmx{config['mem_mb']}M -Djava.io.tmpdir={{output.tempdir}} -Dsnappy.disable=true -Dsamjdk.snappy.disable=true'"
						" FastqToSam"
						" -F1 {input.reads_1}"
						" -F2 {input.reads_2}"
						f" --SAMPLE_NAME {config['alignment_id']}"
						" -O {output.alignments}"
						" --TMP_DIR {output.tempdir}"
						" --USE_JDK_DEFLATER true"
						" --USE_JDK_INFLATER true"


rule merge_aligned_unaligned:
	message:			"Merging aligned and unaligned reads"
	conda:				"../environments/gatk.yaml"
	input:
		aligned			= "alignments/{alignments}.bam",
		unaligned		= f"gatk/{config['alignment_id']}.unaligned.bam",
		reference		= config['reference']
	output:
		alignments		= "gatk/{alignments}.merged.bam",
		tempdir			= temp(directory("temp/gatk/merge_aligned_unaligned/{alignments}"))
	benchmark:			"benchmark/gatk_merge_aligned_unaligned/{alignments}"
	shell:				"mkdir -p {output.tempdir} && gatk"
						f" --java-options '-Xmx{config['mem_mb']}M -Djava.io.tmpdir={{output.tempdir}} -Dsnappy.disable=true -Dsamjdk.snappy.disable=true'"
						" MergeBamAlignment"
						" --TMP_DIR {output.tempdir}"
						" --ATTRIBUTES_TO_RETAIN X0"
						" --ALIGNED_BAM {input.aligned}"
						" --UNALIGNED_BAM {input.unaligned}"
						" --OUTPUT {output.alignments}"
						" --REFERENCE_SEQUENCE {input.reference}"
						" --PAIRED_RUN true"
						" --SORT_ORDER queryname"
						" --IS_BISULFITE_SEQUENCE false"
						" --ALIGNED_READS_ONLY false"
						" --CLIP_ADAPTERS false"
						" --ADD_MATE_CIGAR true"
						" --MAX_INSERTIONS_OR_DELETIONS -1"
						" --PRIMARY_ALIGNMENT_STRATEGY MostDistant"
						" --UNMAPPED_READ_STRATEGY COPY_TO_TAG"
						" --ALIGNER_PROPER_PAIR_FLAGS true"
						" --UNMAP_CONTAMINANT_READS true"
						" --USE_JDK_DEFLATER true"
						" --USE_JDK_INFLATER true"


rule deduplicate:
	message:			"Deduplicating the alignments"
	conda:				"../environments/gatk.yaml"
	input:				"{alignments}.merged.bam"
	output:
		alignments		= "{alignments}.deduplicated.bam",
		dedup_metrics	= "{alignments}.deduplication-metrics",
		tempdir			= temp(directory("temp/gatk/deduplicate/{alignments}"))
	benchmark:			"benchmark/gatk_deduplicate/{alignments}"
	shell:				"mkdir -p {output.tempdir} && gatk"
						f" --java-options '-Xmx{config['mem_mb']}M -Djava.io.tmpdir={{output.tempdir}} -Dsnappy.disable=true -Dsamjdk.snappy.disable=true'"
						" MarkDuplicates"
						" --TMP_DIR {output.tempdir}"
						" --INPUT {input}"
						" --OUTPUT {output.alignments}"
						" --METRICS_FILE {output.dedup_metrics}"
						" --VALIDATION_STRINGENCY SILENT"
						" --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500"
						" --CREATE_MD5_FILE true"
						" --USE_JDK_DEFLATER true"
						" --USE_JDK_INFLATER true"


rule set_nm_uq:
	message:			"Setting the NM and UQ tags"
	conda:				"../environments/gatk.yaml"
	input:
		alignments		= "{alignments}.deduplicated.bam",
		reference		= config['reference']
	output:
		alignments		= "{alignments}.nm-uq-set.bam",
		tempdir			= temp(directory("temp/gatk/set_nm_uq/{alignments}"))
	benchmark:			"benchmark/gatk_set_nm_uq/{alignments}"
	shell:				"mkdir -p {output.tempdir} && gatk"
						f" --java-options '-Xmx{config['mem_mb']}M -Djava.io.tmpdir={{output.tempdir}} -Dsnappy.disable=true -Dsamjdk.snappy.disable=true'"
						" SetNmAndUqTags"
						" --TMP_DIR {output.tempdir}"
						" --INPUT {input.alignments}"
						" --OUTPUT {output.alignments}"
						" --CREATE_INDEX true"
						" --CREATE_MD5_FILE true"
						" --REFERENCE_SEQUENCE {input.reference}"
						" --USE_JDK_DEFLATER true"
						" --USE_JDK_INFLATER true"


rule call:
	message:			"Calling variants with GATK"
	conda:				"../environments/gatk.yaml"
	input:
		alignments		= f"gatk/{config['alignment_id']}.{{wf}}.nm-uq-set.gatk-sorted-coordinate.bam",
		reference		= config['reference'],
		faidx			= f"{config['reference']}.fai",
		ref_dict		= f"{config['reference'].rstrip('.fa')}.dict",
		regions			= "regions/{regions}.bed"
	output:
		variants		= f"variants/{config['alignment_id']}.{{wf}}.gatk.{{regions}}.vcf.gz",
		tempdir			= temp(directory(f"temp/gatk/call/{config['alignment_id']}.{{wf}}.{{regions}}"))
	benchmark:			f"benchmark/gatk_call/{config['alignment_id']}.{{wf}}.{{regions}}"
	shell:				"mkdir -p {output.tempdir} && gatk"
						f" --java-options '-Xmx{config['mem_mb']}M -Djava.io.tmpdir={{output.tempdir}} -Dsnappy.disable=true -Dsamjdk.snappy.disable=true'"
						" HaplotypeCaller"
						f" -ploidy {config['ploidy']}"
						" -R {input.reference}"
						" -I {input.alignments}"
						" -O {output.variants}"
						" -pairHMM LOGLESS_CACHING"
						" -L {input.regions}"
						" --use-jdk-deflater"
						" --use-jdk-inflater"
