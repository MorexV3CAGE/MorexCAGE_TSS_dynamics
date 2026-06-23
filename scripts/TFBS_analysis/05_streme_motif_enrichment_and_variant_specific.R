# ============================================================
# 05) STREME motif enrichment and variant-specific analysis
#
# This consolidated script:
#   1) maps FIMO hits for selected STREME motifs to FASTA windows,
#   2) tests motif enrichment in positive vs background windows,
#   3) tests motif presence in long vs short windows,
#   4) tests whether motifs are variant-specific or shared between
#      long and short promoter variants of the same gene.
#
# Run after:
#   bash 04_scan_streme_motifs_with_fimo.sh
#
# Run:
#   Rscript 05_streme_motif_enrichment_and_variant_specific.R
# ============================================================



# ============================================================
# 0) Libraries
# ============================================================

library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(purrr)

if (requireNamespace("conflicted", quietly = TRUE)) {
  conflicted::conflicts_prefer(base::intersect)
  conflicted::conflicts_prefer(base::setdiff)
  conflicted::conflicts_prefer(base::union)
  conflicted::conflicts_prefer(dplyr::filter)
  conflicted::conflicts_prefer(dplyr::select)
  conflicted::conflicts_prefer(dplyr::mutate)
  conflicted::conflicts_prefer(dplyr::summarise)
  conflicted::conflicts_prefer(dplyr::arrange)
  conflicted::conflicts_prefer(dplyr::count)
  conflicted::conflicts_prefer(dplyr::distinct)
  conflicted::conflicts_prefer(dplyr::first)
}


# ============================================================
# 1) Settings
# ============================================================

min_any_support_for_reporting <- 10L


# ============================================================
# 2) Paths
# ============================================================

project_dir <- Sys.getenv("PROJECT_DIR", unset = ".")
results_dir <- Sys.getenv("RESULTS_DIR", unset = file.path(project_dir, "results"))

streme_discovery_dir <- Sys.getenv(
  "STREME_DISCOVERY_OUTPUT_DIR",
  unset = file.path(results_dir, "streme_motif_discovery")
)

input_dir <- Sys.getenv(
  "INPUT_DIR",
  unset = file.path(streme_discovery_dir, "fasta")
)

pos_fasta <- Sys.getenv(
  "POSITIVE_FASTA",
  unset = file.path(input_dir, "positive_shift_promoter_windows_50_150bp_up.fa")
)

bg_fasta <- Sys.getenv(
  "BACKGROUND_FASTA",
  unset = file.path(input_dir, "background_all_cage_promoter_windows_50_150bp_up.fa")
)

out_root <- Sys.getenv(
  "STREME_MOTIF_OUTPUT_DIR",
  unset = file.path(results_dir, "streme_motif_analysis")
)

table_dir <- file.path(out_root, "tables")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

mapped_fimo_file <- file.path(
  table_dir,
  "01_fimo_streme_motif_hits_mapped_to_windows.tsv"
)

fimo_pos_file <- file.path(
  out_root,
  "fimo",
  "fimo_streme_motifs_positive",
  "fimo.tsv"
)

fimo_bg_file <- file.path(
  out_root,
  "fimo",
  "fimo_streme_motifs_background",
  "fimo.tsv"
)

for (f in c(pos_fasta, bg_fasta)) {
  if (!file.exists(f)) stop("Missing required FASTA file: ", f)
}

if (!file.exists(mapped_fimo_file)) {
  for (f in c(fimo_pos_file, fimo_bg_file)) {
    if (!file.exists(f)) stop("Missing required FIMO file: ", f)
  }
}

message("Using positive FASTA: ", pos_fasta)
message("Using background FASTA: ", bg_fasta)
message("Using STREME motif output directory: ", out_root)


# ============================================================
# 3) Helper functions
# ============================================================

read_fasta_headers <- function(path, set_name) {
  hdr <- readLines(path, warn = FALSE)
  hdr <- hdr[startsWith(hdr, ">")]
  hdr <- sub("^>", "", hdr)

  tibble(set_name = set_name, sequence_name_header = hdr) %>%
    mutate(
      fasta_name = str_replace(.data$sequence_name_header, "::.*$", ""),
      coord_part = if_else(
        str_detect(.data$sequence_name_header, "::"),
        str_replace(.data$sequence_name_header, "^.*::", ""),
        .data$sequence_name_header
      ),
      seqnames = str_match(.data$coord_part, "^([^:]+):")[, 2],
      fasta_start0 = suppressWarnings(as.integer(str_match(.data$coord_part, ":(\\d+)-")[, 2])),
      fasta_end0 = suppressWarnings(as.integer(str_match(.data$coord_part, "-(\\d+)")[, 2])),
      fasta_strand = str_match(.data$coord_part, "\\(([+-])\\)")[, 2],
      geneId = str_match(.data$fasta_name, "^[^|]+\\|([^|]+)\\|")[, 2],
      side = str_match(.data$fasta_name, "^[^|]+\\|[^|]+\\|([^|]+)\\|")[, 2],
      tss = suppressWarnings(as.integer(str_match(.data$fasta_name, "TSS(\\d+)")[, 2])),
      fasta_index = row_number()
    )
}

read_fimo <- function(path, set_name) {
  x <- readr::read_tsv(path, comment = "#", show_col_types = FALSE, progress = FALSE)
  if (nrow(x) == 0) return(tibble())

  names(x) <- make.names(names(x))

  x %>%
    mutate(
      set_name = set_name,
      fimo_hit_id = row_number(),
      sequence_name = as.character(.data$sequence_name),
      fasta_name_from_fimo = str_replace(.data$sequence_name, "::.*$", ""),
      motif_id = as.character(.data$motif_id),
      motif_rank = suppressWarnings(as.integer(str_match(.data$motif_id, "^(\\d+)-")[, 2])),
      motif_consensus = str_replace(.data$motif_id, "^\\d+-", ""),
      fimo_start = as.integer(.data$start),
      fimo_stop = as.integer(.data$stop),
      fimo_strand = as.character(.data$strand),
      fimo_score = as.numeric(.data$score),
      fimo_pvalue = as.numeric(.data$p.value),
      fimo_qvalue = as.numeric(.data$q.value),
      matched_sequence = as.character(.data$matched_sequence)
    )
}

safe_fisher_p <- function(a, b, c, d, alternative = "greater") {
  tryCatch(
    fisher.test(matrix(c(a, b, c, d), nrow = 2, byrow = TRUE), alternative = alternative)$p.value,
    error = function(e) NA_real_
  )
}

safe_fisher_or <- function(a, b, c, d, alternative = "greater") {
  tryCatch(
    unname(fisher.test(matrix(c(a, b, c, d), nrow = 2, byrow = TRUE), alternative = alternative)$estimate),
    error = function(e) NA_real_
  )
}

safe_mcnemar_p <- function(long_only, short_only) {
  if (is.na(long_only) || is.na(short_only)) return(NA_real_)
  if ((long_only + short_only) == 0) return(1)
  tryCatch(
    stats::mcnemar.test(
      matrix(c(0, long_only, short_only, 0), nrow = 2, byrow = TRUE),
      correct = TRUE
    )$p.value,
    error = function(e) NA_real_
  )
}

hyper_p_shared_enriched <- function(shared, n_long, n_short, n_total) {
  stats::phyper(
    q = shared - 1,
    m = n_long,
    n = n_total - n_long,
    k = n_short,
    lower.tail = FALSE
  )
}

hyper_p_variant_specific_enriched <- function(shared, n_long, n_short, n_total) {
  stats::phyper(
    q = shared,
    m = n_long,
    n = n_total - n_long,
    k = n_short,
    lower.tail = TRUE
  )
}

hyper_expected_shared <- function(n_long, n_short, n_total) {
  n_long * n_short / n_total
}

hyper_sd_shared <- function(n_long, n_short, n_total) {
  if (n_total <= 1) return(NA_real_)
  k <- n_short
  K <- n_long
  N <- n_total
  sqrt(k * (K / N) * (1 - K / N) * ((N - k) / (N - 1)))
}

map_fimo_to_windows <- function(fimo_df, header_df) {
  if (nrow(fimo_df) == 0) return(fimo_df)

  fimo_df <- fimo_df %>%
    mutate(
      hit_start0_modeA = .data$fimo_start - 1L,
      hit_end0_modeA = .data$fimo_stop,
      hit_start0_modeB = .data$fimo_start,
      hit_end0_modeB = .data$fimo_stop
    )

  mapped <- vector("list", nrow(fimo_df))

  for (i in seq_len(nrow(fimo_df))) {
    hit <- fimo_df[i, ]

    direct <- header_df %>%
      filter(
        .data$sequence_name_header == hit$sequence_name |
          .data$fasta_name == hit$fasta_name_from_fimo
      )

    if (nrow(direct) > 0) {
      chosen <- direct[1, ]
      mapped[[i]] <- bind_cols(
        hit,
        chosen %>%
          select(
            header_sequence_name = .data$sequence_name_header,
            fasta_name = .data$fasta_name,
            seqnames,
            fasta_start0,
            fasta_end0,
            fasta_strand,
            geneId,
            side,
            tss,
            fasta_index
          )
      ) %>%
        mutate(
          join_method = "direct_header_or_fasta_name",
          hit_start0 = .data$fasta_start0 + .data$fimo_start - 1L,
          hit_end0 = .data$fasta_start0 + .data$fimo_stop
        )
      next
    }

    chr <- hit$sequence_name
    candidates <- header_df %>% filter(.data$seqnames == chr)

    candA <- candidates %>%
      filter(
        .data$fasta_start0 <= hit$hit_start0_modeA,
        .data$fasta_end0 >= hit$hit_end0_modeA
      )

    if (nrow(candA) > 0) {
      chosen <- candA[1, ]
      mapped[[i]] <- bind_cols(
        hit,
        chosen %>%
          select(
            header_sequence_name = .data$sequence_name_header,
            fasta_name = .data$fasta_name,
            seqnames,
            fasta_start0,
            fasta_end0,
            fasta_strand,
            geneId,
            side,
            tss,
            fasta_index
          )
      ) %>%
        mutate(
          join_method = "genomic_coordinate_modeA_1based_to_0based",
          hit_start0 = .data$hit_start0_modeA,
          hit_end0 = .data$hit_end0_modeA
        )
      next
    }

    candB <- candidates %>%
      filter(
        .data$fasta_start0 <= hit$hit_start0_modeB,
        .data$fasta_end0 >= hit$hit_end0_modeB
      )

    if (nrow(candB) > 0) {
      chosen <- candB[1, ]
      mapped[[i]] <- bind_cols(
        hit,
        chosen %>%
          select(
            header_sequence_name = .data$sequence_name_header,
            fasta_name = .data$fasta_name,
            seqnames,
            fasta_start0,
            fasta_end0,
            fasta_strand,
            geneId,
            side,
            tss,
            fasta_index
          )
      ) %>%
        mutate(
          join_method = "genomic_coordinate_modeB_as_0based",
          hit_start0 = .data$hit_start0_modeB,
          hit_end0 = .data$hit_end0_modeB
        )
      next
    }

    mapped[[i]] <- hit %>%
      mutate(
        header_sequence_name = NA_character_,
        fasta_name = NA_character_,
        seqnames = NA_character_,
        fasta_start0 = NA_integer_,
        fasta_end0 = NA_integer_,
        fasta_strand = NA_character_,
        geneId = NA_character_,
        side = NA_character_,
        tss = NA_integer_,
        fasta_index = NA_integer_,
        join_method = NA_character_,
        hit_start0 = NA_integer_,
        hit_end0 = NA_integer_
      )
  }

  bind_rows(mapped)
}


# ============================================================
# 4) Read FASTA headers
# ============================================================

pos_headers <- read_fasta_headers(pos_fasta, "positive")
bg_headers <- read_fasta_headers(bg_fasta, "background")

pos_total <- nrow(pos_headers)
bg_total <- nrow(bg_headers)


# ============================================================
# 5) Load or build mapped FIMO table
# ============================================================

if (file.exists(mapped_fimo_file)) {
  message("Loading existing mapped FIMO table: ", mapped_fimo_file)
  fimo_all <- readr::read_tsv(mapped_fimo_file, show_col_types = FALSE)
} else {
  message("Mapped FIMO table not found; mapping FIMO hits to windows now.")

  fimo_pos <- read_fimo(fimo_pos_file, "positive")
  fimo_bg <- read_fimo(fimo_bg_file, "background")

  fimo_pos_annot <- map_fimo_to_windows(fimo_pos, pos_headers)
  fimo_bg_annot <- map_fimo_to_windows(fimo_bg, bg_headers)

  fimo_all <- bind_rows(fimo_pos_annot, fimo_bg_annot)

  readr::write_tsv(
    fimo_all,
    mapped_fimo_file
  )
}

# Ensure motif_rank/motif_consensus exist if loading an older table.
if (!"motif_rank" %in% colnames(fimo_all)) {
  fimo_all <- fimo_all %>%
    mutate(motif_rank = suppressWarnings(as.integer(str_match(.data$motif_id, "^(\\d+)-")[, 2])))
}
if (!"motif_consensus" %in% colnames(fimo_all)) {
  fimo_all <- fimo_all %>%
    mutate(motif_consensus = str_replace(.data$motif_id, "^\\d+-", ""))
}

fimo_pos_annot <- fimo_all %>% filter(.data$set_name == "positive")
fimo_bg_annot <- fimo_all %>% filter(.data$set_name == "background")

join_diagnostic <- fimo_all %>%
  count(.data$set_name, .data$join_method, is_unmapped = is.na(.data$fasta_name), name = "n_hits") %>%
  arrange(.data$set_name, .data$join_method)

readr::write_tsv(
  join_diagnostic,
  file.path(table_dir, "02_fimo_join_diagnostic.tsv")
)


# ============================================================
# 6) Motif info
# ============================================================

motif_info <- fimo_all %>%
  distinct(.data$motif_id, .data$motif_rank, .data$motif_consensus) %>%
  arrange(.data$motif_rank, .data$motif_id)


# ============================================================
# 7) Positive vs background Fisher per motif
# ============================================================

pos_counts <- fimo_pos_annot %>%
  filter(!is.na(.data$fasta_name)) %>%
  distinct(.data$motif_id, .data$fasta_name) %>%
  count(.data$motif_id, name = "positive_with_motif")

bg_counts <- fimo_bg_annot %>%
  filter(!is.na(.data$fasta_name)) %>%
  distinct(.data$motif_id, .data$fasta_name) %>%
  count(.data$motif_id, name = "background_with_motif")

raw_counts <- fimo_all %>%
  count(.data$motif_id, .data$set_name, name = "raw_fimo_hits") %>%
  pivot_wider(
    names_from = set_name,
    values_from = raw_fimo_hits,
    names_prefix = "raw_hits_",
    values_fill = 0
  )

pos_bg_results <- motif_info %>%
  left_join(pos_counts, by = "motif_id") %>%
  left_join(bg_counts, by = "motif_id") %>%
  left_join(raw_counts, by = "motif_id") %>%
  mutate(
    positive_with_motif = replace_na(.data$positive_with_motif, 0L),
    background_with_motif = replace_na(.data$background_with_motif, 0L),
    raw_hits_positive = replace_na(.data$raw_hits_positive, 0L),
    raw_hits_background = replace_na(.data$raw_hits_background, 0L),
    positive_total = pos_total,
    background_total = bg_total,
    positive_without_motif = .data$positive_total - .data$positive_with_motif,
    background_without_motif = .data$background_total - .data$background_with_motif,
    positive_fraction = .data$positive_with_motif / .data$positive_total,
    background_fraction = .data$background_with_motif / .data$background_total,
    fold_positive_vs_background = .data$positive_fraction / .data$background_fraction,
    fisher_p_positive_enriched = purrr::pmap_dbl(
      list(
        .data$positive_with_motif,
        .data$positive_without_motif,
        .data$background_with_motif,
        .data$background_without_motif
      ),
      ~ safe_fisher_p(..1, ..2, ..3, ..4, alternative = "greater")
    ),
    odds_ratio_positive_vs_background = purrr::pmap_dbl(
      list(
        .data$positive_with_motif,
        .data$positive_without_motif,
        .data$background_with_motif,
        .data$background_without_motif
      ),
      ~ safe_fisher_or(..1, ..2, ..3, ..4, alternative = "greater")
    ),
    fisher_padj_BH_positive_enriched = p.adjust(.data$fisher_p_positive_enriched, method = "BH")
  ) %>%
  arrange(.data$motif_rank, .data$motif_id)

readr::write_tsv(
  pos_bg_results,
  file.path(table_dir, "03_motif_enrichment_positive_vs_background.tsv")
)


# ============================================================
# 8) Long vs short Fisher per motif, window-level
# ============================================================

side_denominator <- pos_headers %>%
  count(.data$side, name = "n_side_total")

long_total <- side_denominator$n_side_total[side_denominator$side == "long"]
short_total <- side_denominator$n_side_total[side_denominator$side == "short"]

side_counts <- fimo_pos_annot %>%
  filter(!is.na(.data$fasta_name), !is.na(.data$side)) %>%
  distinct(.data$motif_id, .data$side, .data$fasta_name) %>%
  count(.data$motif_id, .data$side, name = "n_windows_with_motif") %>%
  pivot_wider(
    names_from = side,
    values_from = n_windows_with_motif,
    values_fill = 0
  )

raw_side_counts <- fimo_pos_annot %>%
  filter(!is.na(.data$side)) %>%
  count(.data$motif_id, .data$side, name = "raw_hits") %>%
  pivot_wider(
    names_from = side,
    values_from = raw_hits,
    names_prefix = "raw_hits_",
    values_fill = 0
  )

side_results <- motif_info %>%
  left_join(side_counts, by = "motif_id") %>%
  left_join(raw_side_counts, by = "motif_id") %>%
  mutate(
    long = replace_na(.data$long, 0L),
    short = replace_na(.data$short, 0L),
    raw_hits_long = replace_na(.data$raw_hits_long, 0L),
    raw_hits_short = replace_na(.data$raw_hits_short, 0L),
    long_total = long_total,
    short_total = short_total,
    long_without_motif = .data$long_total - .data$long,
    short_without_motif = .data$short_total - .data$short,
    long_fraction = .data$long / .data$long_total,
    short_fraction = .data$short / .data$short_total,
    fold_short_vs_long = .data$short_fraction / .data$long_fraction,
    fold_long_vs_short = .data$long_fraction / .data$short_fraction,

    fisher_p_short_enriched_vs_long = purrr::pmap_dbl(
      list(.data$short, .data$short_without_motif, .data$long, .data$long_without_motif),
      ~ safe_fisher_p(..1, ..2, ..3, ..4, alternative = "greater")
    ),
    fisher_p_long_enriched_vs_short = purrr::pmap_dbl(
      list(.data$long, .data$long_without_motif, .data$short, .data$short_without_motif),
      ~ safe_fisher_p(..1, ..2, ..3, ..4, alternative = "greater")
    ),
    fisher_padj_BH_short_enriched_vs_long = p.adjust(.data$fisher_p_short_enriched_vs_long, method = "BH"),
    fisher_padj_BH_long_enriched_vs_short = p.adjust(.data$fisher_p_long_enriched_vs_short, method = "BH")
  ) %>%
  arrange(.data$motif_rank, .data$motif_id)

readr::write_tsv(
  side_results,
  file.path(table_dir, "04_motif_long_vs_short_window_level.tsv")
)


# ============================================================
# 9) Gene set: include genes with >=1 long and >=1 short
# ============================================================

gene_side_window_counts <- pos_headers %>%
  filter(.data$side %in% c("long", "short"), !is.na(.data$geneId)) %>%
  distinct(.data$geneId, .data$side, .data$fasta_name) %>%
  count(.data$geneId, .data$side, name = "n_windows_for_side") %>%
  pivot_wider(
    names_from = side,
    values_from = n_windows_for_side,
    values_fill = 0,
    names_prefix = "n_"
  ) %>%
  mutate(
    n_long = replace_na(.data$n_long, 0L),
    n_short = replace_na(.data$n_short, 0L),
    gene_pair_class = case_when(
      .data$n_long == 1L & .data$n_short == 1L ~ "exactly_one_long_one_short",
      .data$n_long >= 1L & .data$n_short >= 1L ~ "multiple_long_or_short_included",
      .data$n_long >= 1L & .data$n_short == 0L ~ "long_only_gene_no_short_window",
      .data$n_long == 0L & .data$n_short >= 1L ~ "short_only_gene_no_long_window",
      TRUE ~ "no_long_no_short"
    )
  )

long_short_gene_set <- gene_side_window_counts %>%
  filter(.data$n_long >= 1L, .data$n_short >= 1L)

gene_set_summary <- tibble(
  n_positive_windows_total = nrow(pos_headers),
  n_long_windows_total = sum(pos_headers$side == "long", na.rm = TRUE),
  n_short_windows_total = sum(pos_headers$side == "short", na.rm = TRUE),
  n_unique_genes_in_positive_headers = n_distinct(pos_headers$geneId[!is.na(pos_headers$geneId)]),
  n_genes_with_at_least_one_long_and_one_short_window = nrow(long_short_gene_set),
  n_genes_exactly_one_long_one_short = sum(gene_side_window_counts$gene_pair_class == "exactly_one_long_one_short"),
  n_genes_multiple_long_or_short_included = sum(gene_side_window_counts$gene_pair_class == "multiple_long_or_short_included"),
  n_genes_with_long_only_no_short_window = sum(gene_side_window_counts$gene_pair_class == "long_only_gene_no_short_window"),
  n_genes_with_short_only_no_long_window = sum(gene_side_window_counts$gene_pair_class == "short_only_gene_no_long_window")
)

readr::write_tsv(
  gene_side_window_counts,
  file.path(table_dir, "05_gene_long_short_window_counts.tsv")
)

readr::write_tsv(
  gene_set_summary,
  file.path(table_dir, "06_gene_set_summary.tsv")
)


# ============================================================
# 10) Gene-level motif presence: any long vs any short
# ============================================================

motif_gene_side_presence <- fimo_pos_annot %>%
  filter(!is.na(.data$fasta_name), !is.na(.data$geneId), .data$side %in% c("long", "short")) %>%
  distinct(.data$motif_id, .data$geneId, .data$side, .data$fasta_name) %>%
  count(.data$motif_id, .data$geneId, .data$side, name = "n_windows_with_motif_on_side") %>%
  pivot_wider(
    names_from = side,
    values_from = n_windows_with_motif_on_side,
    values_fill = 0,
    names_prefix = "n_motif_windows_"
  ) %>%
  mutate(
    n_motif_windows_long = replace_na(.data$n_motif_windows_long, 0L),
    n_motif_windows_short = replace_na(.data$n_motif_windows_short, 0L),
    long = .data$n_motif_windows_long > 0,
    short = .data$n_motif_windows_short > 0
  )

motif_gene_grid <- tidyr::crossing(
  motif_info %>% select(.data$motif_id, .data$motif_rank, .data$motif_consensus),
  long_short_gene_set %>% select(.data$geneId, .data$n_long, .data$n_short, .data$gene_pair_class)
)

gene_level_detail <- motif_gene_grid %>%
  left_join(
    motif_gene_side_presence,
    by = c("motif_id", "geneId")
  ) %>%
  mutate(
    n_motif_windows_long = replace_na(.data$n_motif_windows_long, 0L),
    n_motif_windows_short = replace_na(.data$n_motif_windows_short, 0L),
    long = replace_na(.data$long, FALSE),
    short = replace_na(.data$short, FALSE),
    motif_presence_class = case_when(
      .data$long & .data$short ~ "shared_both_long_and_short",
      .data$long & !.data$short ~ "long_only",
      !.data$long & .data$short ~ "short_only",
      TRUE ~ "neither"
    ),
    variant_specific = xor(.data$long, .data$short),
    shared = .data$long & .data$short,
    any_variant = .data$long | .data$short
  ) %>%
  arrange(.data$motif_rank, .data$geneId)

readr::write_tsv(
  gene_level_detail,
  file.path(table_dir, "07_motif_presence_by_gene_detail.tsv")
)


# ============================================================
# 11) Gene-level variant-specific vs shared summary + tests
# ============================================================

variant_summary <- gene_level_detail %>%
  group_by(.data$motif_id, .data$motif_rank, .data$motif_consensus) %>%
  summarise(
    gene_total_with_at_least_one_long_and_short = n_distinct(.data$geneId),

    shared_both_long_and_short = sum(.data$shared, na.rm = TRUE),
    long_only = sum(.data$long & !.data$short, na.rm = TRUE),
    short_only = sum(!.data$long & .data$short, na.rm = TRUE),
    variant_specific = sum(.data$variant_specific, na.rm = TRUE),
    neither = sum(!.data$long & !.data$short, na.rm = TRUE),

    genes_with_motif_in_long = sum(.data$long, na.rm = TRUE),
    genes_with_motif_in_short = sum(.data$short, na.rm = TRUE),
    genes_with_motif_in_any_variant = sum(.data$any_variant, na.rm = TRUE),

    fraction_shared_among_all_long_short_genes =
      .data$shared_both_long_and_short / .data$gene_total_with_at_least_one_long_and_short,
    fraction_variant_specific_among_all_long_short_genes =
      .data$variant_specific / .data$gene_total_with_at_least_one_long_and_short,
    fraction_neither_among_all_long_short_genes =
      .data$neither / .data$gene_total_with_at_least_one_long_and_short,

    fraction_shared_among_genes_with_motif_any_variant =
      if_else(.data$genes_with_motif_in_any_variant > 0,
              .data$shared_both_long_and_short / .data$genes_with_motif_in_any_variant,
              NA_real_),

    fraction_variant_specific_among_genes_with_motif_any_variant =
      if_else(.data$genes_with_motif_in_any_variant > 0,
              .data$variant_specific / .data$genes_with_motif_in_any_variant,
              NA_real_),

    expected_shared_under_independent_marginals =
      hyper_expected_shared(
        .data$genes_with_motif_in_long,
        .data$genes_with_motif_in_short,
        .data$gene_total_with_at_least_one_long_and_short
      ),

    sd_shared_under_independent_marginals =
      hyper_sd_shared(
        .data$genes_with_motif_in_long,
        .data$genes_with_motif_in_short,
        .data$gene_total_with_at_least_one_long_and_short
      ),

    shared_observed_minus_expected =
      .data$shared_both_long_and_short - .data$expected_shared_under_independent_marginals,

    shared_observed_over_expected =
      if_else(.data$expected_shared_under_independent_marginals > 0,
              .data$shared_both_long_and_short / .data$expected_shared_under_independent_marginals,
              NA_real_),

    z_shared_vs_independent_marginals =
      if_else(.data$sd_shared_under_independent_marginals > 0,
              .data$shared_observed_minus_expected / .data$sd_shared_under_independent_marginals,
              NA_real_),

    p_shared_enriched_hypergeom =
      hyper_p_shared_enriched(
        .data$shared_both_long_and_short,
        .data$genes_with_motif_in_long,
        .data$genes_with_motif_in_short,
        .data$gene_total_with_at_least_one_long_and_short
      ),

    p_variant_specific_enriched_hypergeom =
      hyper_p_variant_specific_enriched(
        .data$shared_both_long_and_short,
        .data$genes_with_motif_in_long,
        .data$genes_with_motif_in_short,
        .data$gene_total_with_at_least_one_long_and_short
      ),

    binom_p_variant_specific_gt_shared =
      if_else(
        .data$genes_with_motif_in_any_variant > 0,
        purrr::map2_dbl(
          .data$variant_specific,
          .data$genes_with_motif_in_any_variant,
          ~ stats::binom.test(.x, .y, p = 0.5, alternative = "greater")$p.value
        ),
        NA_real_
      ),

    mcnemar_p_long_vs_short_asymmetry =
      safe_mcnemar_p(.data$long_only, .data$short_only),

    passes_min_any_support =
      .data$genes_with_motif_in_any_variant >= min_any_support_for_reporting,

    .groups = "drop"
  ) %>%
  mutate(
    padj_BH_shared_enriched_hypergeom =
      p.adjust(.data$p_shared_enriched_hypergeom, method = "BH"),

    padj_BH_variant_specific_enriched_hypergeom =
      p.adjust(.data$p_variant_specific_enriched_hypergeom, method = "BH"),

    binom_padj_BH_variant_specific_gt_shared =
      p.adjust(.data$binom_p_variant_specific_gt_shared, method = "BH"),

    padj_BH_mcnemar_long_vs_short_asymmetry =
      p.adjust(.data$mcnemar_p_long_vs_short_asymmetry, method = "BH"),

    dominant_pattern_marginal_test = case_when(
      .data$padj_BH_shared_enriched_hypergeom < 0.05 &
        .data$shared_observed_minus_expected > 0 ~ "shared_enriched_vs_marginal_expectation",
      .data$padj_BH_variant_specific_enriched_hypergeom < 0.05 &
        .data$shared_observed_minus_expected < 0 ~ "variant_specific_enriched_vs_marginal_expectation",
      TRUE ~ "not_significant_vs_marginal_expectation"
    )
  ) %>%
  arrange(.data$motif_rank, .data$motif_id)

readr::write_tsv(
  variant_summary,
  file.path(table_dir, "08_motif_variant_specific_shared_summary.tsv")
)


# ============================================================
# 12) Compact reporting table
# ============================================================

reporting_table <- variant_summary %>%
  filter(.data$passes_min_any_support) %>%
  select(
    .data$motif_rank,
    .data$motif_id,
    .data$motif_consensus,
    .data$gene_total_with_at_least_one_long_and_short,
    .data$genes_with_motif_in_any_variant,
    .data$shared_both_long_and_short,
    .data$variant_specific,
    .data$long_only,
    .data$short_only,
    .data$neither,
    .data$fraction_shared_among_genes_with_motif_any_variant,
    .data$fraction_variant_specific_among_genes_with_motif_any_variant,
    .data$binom_p_variant_specific_gt_shared,
    .data$binom_padj_BH_variant_specific_gt_shared,
    .data$expected_shared_under_independent_marginals,
    .data$shared_observed_minus_expected,
    .data$shared_observed_over_expected,
    .data$p_shared_enriched_hypergeom,
    .data$padj_BH_shared_enriched_hypergeom,
    .data$p_variant_specific_enriched_hypergeom,
    .data$padj_BH_variant_specific_enriched_hypergeom,
    .data$mcnemar_p_long_vs_short_asymmetry,
    .data$padj_BH_mcnemar_long_vs_short_asymmetry,
    .data$dominant_pattern_marginal_test
  )

readr::write_tsv(
  reporting_table,
  file.path(table_dir, "09_motif_variant_specific_shared_reporting_table.tsv")
)


# ============================================================
# 13) Console summary
# ============================================================

message("\nDONE STREME motif enrichment and variant-specific analysis")
message("Output: ", table_dir)

message("\nGene set summary:")
print(gene_set_summary)

message("\nPositive vs background Fisher:")
print(pos_bg_results)

message("\nWindow-level long vs short Fisher:")
print(side_results)

message("\nGene-level variant-specific vs shared summary:")
print(variant_summary)

message("\nCompact reporting table:")
print(reporting_table)
