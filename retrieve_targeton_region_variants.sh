#! /bin/bash
#BSUB -q normal
#BSUB -G cellularoperations
#BSUB -J 1_1_A
#BSUB -R "select[mem>5000] rusage[mem=5000] span[hosts=1]"
#BSUB -M 5000
#BSUB -n 4
#BSUB -oo /lustre/scratch117/sciops/team302/targeton_design/targeton_design_v021/subset_variants.o
#BSUB -eo /lustre/scratch117/sciops/team302/targeton_design/targeton_design_v021/subset_variants.err

module load /software/modules/common-apps/bcftools/1.9-220

resources_dir="/lustre/scratch117/sciops/team302/resources"
clinvar_ver="clinvar_"
gene=""
chr="chr"
targeton_design_dir=""
regions_minus_primers="${targeton_design_dir}/intermediate_outputs/${gene}_primers.txt"

wget "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/${clinvar_ver}.vcf.gz"
wget "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/${clinvar_ver}.vcf.gz.tbi"




# Update of chromosome names required for clinvar 1 -> chr1
bcftools annotate --rename-chrs "${resources_dir}/update_chroms.txt" \
"${resources_dir}/${clinvar_ver}.vcf.gz" -O z -o "${resources_dir}/${clinvar_ver}_rename.vcf.gz"
# Reindex clinvar
bcftools index "${resources_dir}/${clinvar_ver}_rename.vcf.gz"


# Subset clinvar regions
bcftools view -R "${regions_minus_primers}" \
"${resources_dir}/${clinvar_ver}_rename.vcf.gz" \
-o "${targeton_design_dir}/intermediate_outputs/variant_sources/${clinvar_ver}_${gene}.vcf"

# rsync /warehouse/mave/sge_production_02/variation_data/gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz ${resources_dir}
# Subset gnomad2 regions
bcftools view -R "${regions_minus_primers}" \
"${resources_dir}/gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz" \
-o "${targeton_design_dir}/intermediate_outputs/variant_sources/gnomAD_v2.1.1_${gene}_hg38.vcf"

# Subset gnomad3 regions
bcftools view -R "${regions_minus_primers}" \
"/lustre/scratch117/core/corebio/gnomad/3.1/variants/genomes/gnomad.genomes.v3.1.sites.${chr}.vcf.bgz" \
-o "${targeton_design_dir}/intermediate_outputs/variant_sources/gnomAD_v3.1_${gene}_hg38.vcf"



# Subset supplementary regions
bcftools sort "${targeton_design_dir}/xxx.VCF" \
-O z -o "${targeton_design_dir}/intermediate_outputs/variant_sources/${gene}_additional_vars.vcf.bgz" 
bcftools index "${targeton_design_dir}/intermediate_outputs/variant_sources/${gene}_additional_vars.vcf.bgz"


bcftools view -R "${regions_minus_primers}" \
"${targeton_design_dir}/intermediate_outputs/variant_sources/${gene}_additional_vars.vcf.bgz" \
 -o "${targeton_design_dir}/intermediate_outputs/variant_sources/${gene}_subset_vars.vcf"

