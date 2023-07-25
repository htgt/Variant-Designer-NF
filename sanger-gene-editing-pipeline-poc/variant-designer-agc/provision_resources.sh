#!/usr/bin/env bash
set -x

CLINVAR_VERSION=clinvar
GNOMAD_CHR=16
MANE_VER=1.0
S3_URI=$1

# remove trailing /
S3_URI=`echo $S3_URI | sed 's:/$::'`

# create temp dir
TMPDIR=$(mktemp -d)

# inputs
# ======
aws s3 cp inputs/CTCF_target_regions_manifest.txt ${S3_URI}/project/GeneEditingPipeline/inputs/VariantDesigner/
aws s3 cp inputs/design_choices_51-63.tsv ${S3_URI}/project/GeneEditingPipeline/inputs/VariantDesigner/

# R
# =
aws s3 cp bin/validation_utils.R ${S3_URI}/project/GeneEditingPipeline/resources/VariantDesigner/R/

# other
# =====
aws s3 cp resources/update_chroms.txt ${S3_URI}/project/GeneEditingPipeline/resources/VariantDesigner/other/


# MANE Summary
# =====
wget --directory-prefix=${TMPDIR}  https://ftp.funet.fi/pub/mirrors/ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_${MANE_VER}/MANE.GRCh38.v${MANE_VER}.summary.txt.gz
aws s3 cp ${TMPDIR}/MANE.GRCh38.v${MANE_VER}.summary.txt.gz ${S3_URI}/project/GeneEditingPipeline/resources/VariantDesigner/mane/

# Assembly file
# =====
wget -q -O ${TMPDIR}/Homo_sapiens_assembly38_chr16.fa "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=CM000678.2&rettype=fasta"
aws s3 cp ${TMPDIR}/Homo_sapiens_assembly38_chr16.fa ${S3_URI}/project/GeneEditingPipeline/resources/VariantDesigner/assembly/

# ClinVar
# =======
wget --directory-prefix=${TMPDIR} https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/${CLINVAR_VERSION}.vcf.gz
wget --directory-prefix=${TMPDIR} https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/${CLINVAR_VERSION}.vcf.gz.tbi

aws s3 cp ${TMPDIR}/${CLINVAR_VERSION}.vcf.gz ${S3_URI}/project/GeneEditingPipeline/resources/VariantDesigner/clinvar/
aws s3 cp ${TMPDIR}/${CLINVAR_VERSION}.vcf.gz.tbi ${S3_URI}/project/GeneEditingPipeline/resources/VariantDesigner/clinvar/

# Gnomad
# ======
aws s3 cp s3://gnomad-public-us-east-1/release/2.1.1/liftover_grch38/vcf/exomes/gnomad.exomes.r2.1.1.sites.${GNOMAD_CHR}.liftover_grch38.vcf.bgz ${S3_URI}/project/GeneEditingPipeline/resources/VariantDesigner/gnomad/exomes/
aws s3 cp s3://gnomad-public-us-east-1/release/2.1.1/liftover_grch38/vcf/exomes/gnomad.exomes.r2.1.1.sites.${GNOMAD_CHR}.liftover_grch38.vcf.bgz.tbi ${S3_URI}/project/GeneEditingPipeline/resources/VariantDesigner/gnomad/exomes/

aws s3 cp s3://gnomad-public-us-east-1/release/3.1/vcf/genomes/gnomad.genomes.v3.1.sites.chr${GNOMAD_CHR}.vcf.bgz ${S3_URI}/project/GeneEditingPipeline/resources/VariantDesigner/gnomad/genomes/
aws s3 cp s3://gnomad-public-us-east-1/release/3.1/vcf/genomes/gnomad.genomes.v3.1.sites.chr${GNOMAD_CHR}.vcf.bgz.tbi ${S3_URI}/project/GeneEditingPipeline/resources/VariantDesigner/gnomad/genomes/


# remove temp dir
rm -r $TMPDIR
