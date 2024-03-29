#!/usr/bin/env Rscript

# get command line arguments
args = commandArgs(trailingOnly=TRUE)
validation_utils_r <- args[1]
glue_file <- args[2]
targets_tsv <- args[3]
targetons_tsv <- args[4]

################################################################################
# Setup depedencies ------------------------------------------------------------
################################################################################
library(validate)
library(tidyverse)
library(glue)
library(here)
options(scipen = 100)
library(Biostrings)
library(vroom)
library(cli)
print(getwd())
print(dir())
source(validation_utils_r)
mane_version <- "1.0"

mane_summary_data <- read_tsv(
  glue(
    glue_file
  ),
  col_types = cols()
)

range <- "63-70"

################################################################################
# Read in data -----------------------------------------------------------------
################################################################################
print("Read in data")
# replace here("") in these:
print("Reading targets...")
targets <- read_tsv(targets_tsv) 
print("done")
print("Reading targetons...")
targetons <- read_tsv(targetons_tsv)
print("done")
# Add a list column to targets detailing split and regular targetons
# targets[["Targeton_ID"]] <- c(
#   list(targetons[["Targeton_ID"]][1:4]),
#   targetons[["Targeton_ID"]][5:13]
# )
print("################################################################################")
print("Head targets:")
head(targets)
print("Head targetons:")
head(targetons)
print("################################################################################")
################################################################################
# Wrangle datasets -------------------------------------------------------------
################################################################################
print("Get strand for each gene...")
# Get strand for each gene
targets <- targets |>
  mutate(ref_strand = get_gene_strand(Gene)) |>
  unnest(c(`Targeton_ID`))
print("done")
print("Autofill LibAmp primer info...")
# Autofill LibAmp primer info
targetons <- targetons |>
  separate("CDS_position (with splits)", 
    sep = "[:-]",
    into = c("ref_chr", "r2_start", "r2_end")
  ) |>
  mutate(
    `Forward Primer Length` = str_length(`Forward Primer`),
    `Reverse Primer Length` = str_length(`Reverse Primer`),
    `Amplicon length` = `Amplicon End` - `Amplicon Start`,
    `ref_start` = `Amplicon Start`,
    `ref_end` = `Amplicon End`,
    `r2_start` = as.integer(`r2_start`),
    `r2_end` = as.integer(`r2_end`)
  )
print("done")
################################################################################
# Merge targets and targetons --------------------------------------------------
################################################################################
print("B")
combined_data <- left_join(targets,
  targetons,
  by = c("Targeton_ID", "Gene") #only design choice has both targeton and gene 
)

# Calculate the R1 and R3 regions for each targeton
data_plus_ext_vec <- combined_data |>
  mutate(
    poss_5_bases = r2_start - `Amplicon Start` - `Forward Primer Length`,
    poss_3_bases = `Amplicon End` - r2_end - `Reverse Primer Length`
  ) |>
  rowwise() |>
  mutate(
    R1 = case_when(
      ref_strand == "+" ~ min(`Maximum upstream (Default:25)`, poss_5_bases),
      ref_strand == "-" ~ min(`Maximum downstream (Default:15)`, poss_5_bases)
    ),
    R3 = case_when(
      ref_strand == "+" ~ min(`Maximum downstream (Default:15)`, poss_3_bases),
      ref_strand == "-" ~ min(`Maximum upstream (Default:25)`, poss_3_bases)
    ),
    ext_vector = case_when(
      is.na(split_pos) ~ paste(R1, R3, sep = ","),
      split_pos == "five" ~ paste(R1, 0, sep = ","),
      split_pos == "three" ~ paste(0, R3, sep = ","),
      split_pos == "internal" ~ paste(0, 0, sep = ",")
    )
  ) |>
  mutate(`Output targeton coordinates` = paste0(
    ref_chr, ":", r2_start - R1,
    "-", r2_end + R3
  ))

# Calculate the action vector that will be applied to R1, R2 and R3
data_plus_act_vec <- data_plus_ext_vec |>
  rowwise() |>
  mutate(action_vector = calculate_action_vector(
    strand = ref_strand,
    R1 = R1,
    us = `Mutation to be include in upstream (default: snv)`,
    ds = `Mutation to be include in downstream (default: snv)`,
    r2_mutator = `Mutation to be include in R2 (default: snvre,inframe)`,
    split_pos = split_pos
  ))

# Pull sgRNA data from WGE and create the sgRNA vector which IDs a lib
data_plus_sgrna_info <- data_plus_act_vec |>
  mutate(
    sgRNA_sequence = retrieve_wge_data(sgRNA_id)[["grna_seq"]],
    sgRNA_start = retrieve_wge_data(sgRNA_id)[["grna_start"]],
    sgRNA_strand = retrieve_wge_data(sgRNA_id)[["grna_strand"]],
    sgRNA_off_targets = retrieve_wge_data(sgRNA_id)[["grna_off_targets"]],
    .after = sgRNA_id
  ) |>
  mutate(sgrna_vector = paste0("sgRNA_", sgRNA_id), .after = action_vector)


data_plus_sgrna_info[["Region that must be  included in this SGE if we need to split the exon"]] <- NA
################################################################################
# Write output files -----------------------------------------------------------
################################################################################
print("C")
output_manifest <- data_plus_sgrna_info |>
  select(
    Gene,
    Exon_position,
    `Coding Exon Size(bp)`,
    `If the coding Exon>240bp, Do you want to cover the whole exon? (Yes/No, include the selected region only)`,
    `Region that must be  included in this SGE if we need to split the exon`,
    # `PAM_VCF`,
    sgRNA_id,
    sgRNA_sequence,
    sgRNA_start,
    sgRNA_strand,
    sgRNA_off_targets,
    # PAM_VCF,
    `Forward Primer`,
    `Forward Primer Length`,
    `Reverse Primer`,
    `Reverse Primer Length`,
    `Amplicon Start`,
    `Amplicon End`,
    ref_chr,
    ref_strand,
    ref_start,
    ref_end,
    r2_start,
    r2_end,
    ext_vector,
    action_vector,
    sgrna_vector,
    `Output targeton coordinates`,
    Targeton_ID
  )

write_tsv(output_manifest, file = here(glue("targetons_{range}_manifest.tsv")))
print("D")

dispr_input <- output_manifest |>
  select(`Forward Primer`, `Reverse Primer`, ref_chr)
print("E")
write_tsv(dispr_input, "x_primers.tsv")
print("F")
