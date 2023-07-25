#!/bin/bash

# targeton_design_dir=""
regions_minus_primers=$1
gene=$2 # ""
chr=$3 # "chr"
update_chroms=$4
query_chr=$5
clinvar_ver=$6

# Step 1
# Update of chromosome names required for clinvar 1 -> chr1
bcftools annotate --rename-chrs "${update_chroms}" \
${clinvar_ver}.vcf.gz -O z \
-o "${clinvar_ver}_rename.vcf.gz" # outfile1
# Reindex clinvar
bcftools index ${clinvar_ver}_rename.vcf.gz # outfile2

# Step 2
# Subset clinvar regions
bcftools view -R ${regions_minus_primers} \
${clinvar_ver}_rename.vcf.gz \
-o ${clinvar_ver}_${gene}.vcf # outfile3

# Step 3
# Subset gnomad2 regions
bcftools view -R ${regions_minus_primers} \
gnomad.exomes.r2.1.1.sites.${query_chr}.liftover_grch38.vcf.bgz \
-o gnomAD_v2.1.1_${gene}_hg38.vcf # outfile4

# Step 4
# Subset gnomad3 regions
bcftools view -R ${regions_minus_primers} \
gnomad.genomes.v3.1.sites.${chr}.vcf.bgz \
-o gnomAD_v3.1_${gene}_hg38.vcf # outfile5

# Step 5
# Subset supplementary regions
# bcftools sort "${targeton_design_dir}/xxx.VCF" \
# -O z -o "${targeton_design_dir}/intermediate_outputs/variant_sources/${gene}_additional_vars.vcf.bgz" 
# bcftools index "${targeton_design_dir}/intermediate_outputs/variant_sources/${gene}_additional_vars.vcf.bgz"

# bcftools view -R "${regions_minus_primers}" \
# "${targeton_design_dir}/intermediate_outputs/variant_sources/${gene}_additional_vars.vcf.bgz" \
#  -o "${targeton_design_dir}/intermediate_outputs/variant_sources/${gene}_subset_vars.vcf"
