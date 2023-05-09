#!/usr/bin/env nextflow

// Define variables from input
chrom_file = params.chrom_file
gene_name = params.gene_name
chromosome = params.chromosome
query_chromosome = chromosome.minus("chr")

clinvar_ver = params.clinvar_ver
resource_dir = params.resource_dir

// Stage input files
targets = Channel.fromPath(params.targets)
targetons = Channel.fromPath(params.targetons)

// Stage input resources
validation_utils = Channel.fromPath(params.validation_utils)
glue_file = Channel.fromPath(params.glue_file)
update_chroms = Channel.fromPath(params.update_chroms)
clinvar = Channel.fromPath("${resource_dir}/clinvar/${clinvar_ver}.vcf.gz")
clinvar_tbi = Channel.fromPath("${resource_dir}/clinvar/${clinvar_ver}.vcf.gz.tbi")
gnomad_v2 = Channel.fromPath("${resource_dir}/gnomad/exomes/gnomad.exomes.r2.1.1.sites.${query_chromosome}.liftover_grch38.vcf.bgz")
gnomad_v2_tbi = Channel.fromPath("${resource_dir}/gnomad/exomes/gnomad.exomes.r2.1.1.sites.${query_chromosome}.liftover_grch38.vcf.bgz.tbi")
gnomad_v3 = Channel.fromPath("${resource_dir}/gnomad/genomes/gnomad.genomes.v3.1.sites.${chromosome}.vcf.bgz")
gnomad_v3_tbi = Channel.fromPath("${resource_dir}/gnomad/genomes/gnomad.genomes.v3.1.sites.${chromosome}.vcf.bgz.tbi")

process AutofillTargeton {
    publishDir "${params.output_dir}"
    container = "${params.container_registry}/gene-editing-pipeline/autofill-targeton:1.0.0"
    
    input:
        path validation_utils
        path glue_file
        path targets
        path targetons

    output:
        path 'x_primers.tsv', emit:primer_file
    
    script:
        """
        autofill_targeton_manifest.R $validation_utils $glue_file $targets $targetons
        """
}


process GenerateBedFiles {
    publishDir "${params.output_dir}"
    container = "${params.container_registry}/gene-editing-pipeline/generate-bedfiles:1.0.0"
    input:
        path(chrom_file)
        path(primer_file)
        val(gene_name)

    output:
        path "${gene_name}_internal.bed", emit:primer_bed
    
    script:
    """
    generate_bedfiles.sh $primer_file $gene_name $chrom_file
    """
}


process RetrieveTargetons {
    publishDir "${params.output_dir}"
    container = "${params.container_registry}/gene-editing-pipeline/retrieve-targetons:1.0.0"

    input:
        path(primer_bed)
        val(gene_name)
        val(chromosome)
        path(update_chroms)
        val(query_chromosome)
        val(clinvar_ver)
        path(clinvar)
        path(clinvar_tbi)
        path(gnomad_v2)
        path(gnomad_v2_tbi)
        path(gnomad_v3)
        path(gnomad_v3_tbi)

    output:
        path "clinvar_rename.vcf.gz"
        path "clinvar_rename.vcf.gz.csi"
        path "clinvar_${gene_name}.vcf"
        path "gnomAD_v2.1.1_${gene_name}_hg38.vcf"
        path "gnomAD_v3.1_${gene_name}_hg38.vcf", emit:gnomad3_out
    
    script:
    """
    retrieve_targeton_region_variants.sh ${primer_bed} ${gene_name} ${chromosome} ${update_chroms} ${query_chromosome} ${clinvar_ver}
    """
}


workflow {
  AutofillTargeton(validation_utils, glue_file, targets, targetons)
  GenerateBedFiles(chrom_file, AutofillTargeton.out.primer_file, gene_name)
  RetrieveTargetons(GenerateBedFiles.out.primer_bed, gene_name, chromosome, update_chroms, query_chromosome, clinvar_ver, clinvar, clinvar_tbi, gnomad_v2, gnomad_v2_tbi, gnomad_v3, gnomad_v3_tbi)
}
