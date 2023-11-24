#!/bin/bash

set -euxo pipefail

# hs37d5
mkdir -p reference
pushd reference
wget -c https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
if [ ! -e hs37d5.fa ]
then
	gunzip -k hs37d5.fa.gz
fi
popd

# Truthset for NA24385
mkdir -p truth
pushd truth
wget -c https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NIST_SV_v0.6/HG002_SVs_Tier1_v0.6.vcf.gz
wget -c https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NIST_SV_v0.6/HG002_SVs_Tier1_v0.6.bed
popd

# Known variants
mkdir -p known-variants
pushd known-variants
wget -c https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
# Unfortunately we currently need the uncompressed variants.
if [ ! -e ALL.chr1.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf ]
then
	gunzip -k ALL.chr1.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
fi
popd

# Reads
mkdir -p reads
pushd reads
wget -c https://s3-us-west-2.amazonaws.com/human-pangenomics/NHGRI_UCSC_panel/HG002/hpp_HG002_NA24385_son_v1/ILMN/downsampled/HG002_HiSeq30x_subsampled_R1.fastq.gz
wget -c https://s3-us-west-2.amazonaws.com/human-pangenomics/NHGRI_UCSC_panel/HG002/hpp_HG002_NA24385_son_v1/ILMN/downsampled/HG002_HiSeq30x_subsampled_R2.fastq.gz
popd
