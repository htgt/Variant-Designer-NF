revcomp <- function(x) {
  as.character(reverseComplement(DNAStringSet(x)))
}

find_chr_length <- function(x) {
  hg38_chrom_lengths |>
    filter(chrom == x) |>
    pull(length)
}

get_mane_data <- function(x) {
  mane_summary_data |> filter(Ensembl_nuc == x)
}

get_gene_strand <- function(x) {
  mane_summary_data |>
    filter(symbol == x) |>
    pull(chr_strand)
}


get_ens_data <- function(mane_data, x) {
  mane_data |>
    mutate(attr_dat = str_split(attribute, pattern = ";")) |>
    unnest(attr_dat) |>
    filter(grepl(attr_dat, pattern = x))
}

calculate_action_vector <- function(strand, R1, us, ds, r2_mutator, split_pos) {
  rep_str <- ifelse(R1 %% 2 == 0, "2del0", "2del1")
  if (strand == "+") {
    r1_mutator <- str_replace(us, pattern = "2del", replacement = rep_str)
    r3_mutator <- str_replace(ds, pattern = "2del", replacement = "2del0")
  } else {
    r1_mutator <- str_replace(ds, pattern = "2del", replacement = rep_str)
    r3_mutator <- str_replace(us, pattern = "2del", replacement = "2del0")
  }
  if (!is.na(split_pos)) {
    if (split_pos == "five") {
      r3_mutator <- ""
    }
    if (split_pos == "three") {
      r1_mutator <- ""
    }
    if (split_pos == "internal") {
      r1_mutator <- ""
      r3_mutator <- ""
    }
  }

  paste0("(", r1_mutator, ")", ",", "(", r2_mutator, ")", ",", "(", r3_mutator, ")")
}





format_mane_data <- function(x) {
  temp <- x |>
    rowwise() |>
    mutate(col_list = list(str_split(str_remove_all(str_split(attribute, "; ", simplify = TRUE),
      pattern = '\"'
    ), pattern = " "))) |>
    unnest(cols = col_list) |>
    unnest_wider(col_list, names_sep = "_") |>
    pivot_wider(names_from = `col_list_1`, values_from = `col_list_2`) |>
    distinct()

  temp |>
    unnest(cols = colnames(temp)[11:21]) |>
    distinct() |>
    select(-c(attribute, attr_dat)) |>
    filter(feature %in% c("exon", "CDS")) |>
    select(feature, gene_name, seqname, exon_id, start, end, exon_number, transcript_id) |>
    distinct() |>
    pivot_wider(names_from = feature, values_from = c(start, end)) |>
    mutate(
      Exon_position = paste0(seqname, ":", start_exon, "-", end_exon),
      CDS_position = paste0(seqname, ":", start_CDS, "-", end_CDS),
      length = 1 + abs(end_CDS - start_CDS)
    ) |>
    select(gene_name, Exon_position, CDS_position, exon_id, transcript_id, exon_number, length) |>
    mutate(exon_number = as.numeric(exon_number))
}


target_region_exists <- function(gene_name, ex_pos, cds_pos, ex_id, tx_id, ex_number, size, data) {
  mergeable <- merge(tibble(
    gene_name = gene_name,
    exon_number = ex_number,
    Exon_position = ex_pos,
    CDS_position = cds_pos,
    exon_id = ex_id,
    transcript_id = tx_id,
    length = size
  ), y = data) |> nrow()
  mergeable == 1
}



retrieve_wge_data <- function(grna_id) {
  if (!is.na(grna_id)) {
    res <- GET(glue("https://wge.stemcell.sanger.ac.uk/api/crispr_by_id?id={grna_id}&id={grna_id}&species=Grch38"))
    seq <- fromJSON(rawToChar(res$content))[[1]][["seq"]]
    start <- fromJSON(rawToChar(res$content))[[1]][["chr_start"]]
    dir <- fromJSON(rawToChar(res$content))[[1]][["pam_right"]]
    off_targets <- fromJSON(rawToChar(res$content))[[1]][["off_target_summary"]]
    if (dir == 1) {
      sequence <- substr(seq, start = 1, stop = 20)
    } else {
      sequence <- revcomp(substr(seq, start = 4, stop = 23))
    }
    strand <- if_else(dir == 1, "+", "-")

    list(grna_seq = sequence, grna_start = as.numeric(start), grna_strand = strand, grna_off_targets = off_targets)
  } else {
    list(grna_seq = NA, grna_start = NA, grna_strand = NA, off_targets = NA)
  }
}

check_valid_mutators <- function(x, permitted_fields) {
  map_depth(
    map(x,
      str_split,
      pattern = ","
    ),
    ~ all(str_trim(.x) %in% permitted_fields),
    .depth = 2
  ) |>
    purrr::flatten() |>
    purrr::flatten_lgl()
}
