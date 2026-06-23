# ============================================================
# 02) Known TFBS-family module enrichment and variant-specific analysis
#
# Goals:
#   1) summarizes known TFBS overlaps,
#   2) merges overlapping same-family TFBS hits into blocks,
#   3) defines local TFBS-family modules with a 10-50 bp gap,
#   4) tests module enrichment in positive vs background windows,
#   5) tests whether enriched modules are variant-specific or shared
#      between long and short promoter variants of the same gene.
#
# Run after:
#   bash 01_prepare_known_tfbs_window_overlaps.sh
#
# Run:
#   Rscript 02_known_tfbs_module_enrichment_and_variant_specific.R
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

min_gap_bp <- 10L
max_gap_bp <- 50L
ngram_lengths <- c(2L, 3L)

# This is only a reporting flag; all tests are written to the full table.
min_positive_support_for_reporting <- 5L

# If TRUE, unknown family hits are dropped before block/module construction.
drop_unknown_family <- TRUE

# BED column containing motif ID/name. Standard BED name column = 4.
tfbs_motif_col <- 4L


# ============================================================
# 2) Locate manifest from script 01
# ============================================================

project_dir <- Sys.getenv("PROJECT_DIR", unset = ".")
results_dir <- Sys.getenv("RESULTS_DIR", unset = file.path(project_dir, "results"))

out_root <- Sys.getenv(
  "KNOWN_TFBS_OUTPUT_DIR",
  unset = file.path(results_dir, "known_tfbs_module_analysis")
)

table_dir <- file.path(out_root, "tables")
manifest_file <- file.path(table_dir, "00_known_tfbs_window_overlap_manifest.tsv")

if (!file.exists(manifest_file)) {
  stop(
    "Cannot find known TFBS overlap manifest. Did you run script 01 first? Tried:
",
    manifest_file
  )
}

message("Using project_dir: ", project_dir)
message("Using out_root: ", out_root)

manifest <- readr::read_tsv(manifest_file, show_col_types = FALSE)

get_manifest <- function(key, default = NA_character_) {
  val <- manifest$value[manifest$parameter == key]
  if (length(val) == 0 || is.na(val[1]) || val[1] == "") default else val[1]
}

resolve_path <- function(path) {
  if (is.na(path) || path == "") return(path)
  if (file.exists(path)) return(path)
  path
}

pos_windows_bed <- resolve_path(get_manifest("POSITIVE_WINDOWS_BED"))
bg_windows_bed <- resolve_path(get_manifest("BACKGROUND_WINDOWS_BED"))
pos_intersect <- resolve_path(get_manifest("POSITIVE_TFBS_OVERLAPS"))
bg_intersect <- resolve_path(get_manifest("BACKGROUND_TFBS_OVERLAPS"))
family_map_file <- resolve_path(get_manifest("TF_FAMILY_MAP"))
pos_fasta <- resolve_path(get_manifest("POSITIVE_FASTA"))
tfbs_ncol <- as.integer(get_manifest("TFBS_NCOL"))

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

message("Resolved input paths:")
message("  pos_windows_bed: ", pos_windows_bed)
message("  bg_windows_bed:  ", bg_windows_bed)
message("  pos_intersect:   ", pos_intersect)
message("  bg_intersect:    ", bg_intersect)
message("  family_map_file: ", family_map_file)
message("  pos_fasta:       ", pos_fasta)

for (f in c(pos_windows_bed, bg_windows_bed, pos_intersect, bg_intersect, pos_fasta)) {
  if (!file.exists(f)) {
    stop(
      "Missing required file after path resolution: ", f,
      "\nCheck that manifest paths point to accessible files."
    )
  }
}

if (is.na(tfbs_ncol) || tfbs_ncol < 4L) {
  stop("Invalid TFBS_NCOL from manifest: ", get_manifest("TFBS_NCOL"))
}

message("TFBS_NCOL: ", tfbs_ncol)
message("family_map_file: ", family_map_file)


# ============================================================
# 3) Helper functions
# ============================================================

read_windows_bed <- function(path) {
  readr::read_tsv(
    path,
    col_names = c("seqnames", "window_start0", "window_end0", "fasta_name", "set_name", "side", "geneId", "tss", "window_strand"),
    show_col_types = FALSE
  ) %>%
    mutate(
      window_start0 = as.integer(.data$window_start0),
      window_end0 = as.integer(.data$window_end0),
      tss = suppressWarnings(as.integer(.data$tss))
    )
}

read_intersect <- function(path, set_name, tfbs_ncol) {
  x <- readr::read_tsv(path, col_names = FALSE, show_col_types = FALSE, progress = FALSE)

  if (nrow(x) == 0) {
    return(tibble())
  }

  window_offset <- tfbs_ncol

  out <- tibble(
    set_name = set_name,
    tfbs_seqnames = as.character(x[[1]]),
    tfbs_start0 = as.integer(x[[2]]),
    tfbs_end0 = as.integer(x[[3]]),
    motif_id_raw = as.character(x[[tfbs_motif_col]]),

    window_seqnames = as.character(x[[window_offset + 1]]),
    window_start0 = as.integer(x[[window_offset + 2]]),
    window_end0 = as.integer(x[[window_offset + 3]]),
    fasta_name = as.character(x[[window_offset + 4]]),
    window_set_name = as.character(x[[window_offset + 5]]),
    side = as.character(x[[window_offset + 6]]),
    geneId = as.character(x[[window_offset + 7]]),
    tss = suppressWarnings(as.integer(x[[window_offset + 8]])),
    window_strand = as.character(x[[window_offset + 9]])
  ) %>%
    mutate(
      # Clip TFBS to window. This keeps module construction inside one 50-150 bp window.
      hit_start0 = pmax(.data$tfbs_start0, .data$window_start0),
      hit_end0 = pmin(.data$tfbs_end0, .data$window_end0),
      hit_width = .data$hit_end0 - .data$hit_start0,
      motif_id = .data$motif_id_raw
    ) %>%
    filter(.data$hit_width > 0)

  out
}

normalize_id <- function(x) {
  x <- as.character(x)
  x <- str_replace(x, "^MOTIF\\s+", "")
  x <- str_replace(x, "\\s.*$", "")
  x <- str_replace(x, "\\|.*$", "")
  x <- str_replace(x, ";.*$", "")
  x <- str_replace(x, "__.*$", "")
  x
}

normalize_family <- function(x) {
  x <- as.character(x)
  case_when(
    is.na(x) | x == "" ~ NA_character_,
    x %in% c("BBR_BPC", "BBR-BPC") ~ "BBR-BPC",
    TRUE ~ x
  )
}

read_plant_tf_map <- function(path) {
  if (is.na(path) || path == "" || !file.exists(path)) {
    warning("PlantTFDB family map not found: ", path)
    return(tibble(
      TF_ID = character(),
      TF_ID_no_iso = character(),
      Gene_ID = character(),
      Gene_ID_no_iso = character(),
      Family = character()
    ))
  }

  x <- readr::read_tsv(path, show_col_types = FALSE)

  # Try to standardize common PlantTFDB table column names.
  cn <- colnames(x)
  if (!"TF_ID" %in% cn) {
    possible <- cn[str_detect(cn, regex("^TF.*ID$|^TF_ID$|^Name$|^ID$", ignore_case = TRUE))]
    if (length(possible) > 0) names(x)[names(x) == possible[1]] <- "TF_ID"
  }
  if (!"Gene_ID" %in% cn) {
    possible <- colnames(x)[str_detect(colnames(x), regex("^Gene.*ID$|^Gene_ID$|^geneId$", ignore_case = TRUE))]
    if (length(possible) > 0) names(x)[names(x) == possible[1]] <- "Gene_ID"
  }
  if (!"Family" %in% colnames(x)) {
    possible <- colnames(x)[str_detect(colnames(x), regex("family", ignore_case = TRUE))]
    if (length(possible) > 0) names(x)[names(x) == possible[1]] <- "Family"
  }

  required <- c("TF_ID", "Gene_ID", "Family")
  missing <- setdiff(required, colnames(x))
  if (length(missing) > 0) {
    stop(
      "Family map does not contain required columns: ",
      paste(missing, collapse = ", "),
      "\nColumns found: ",
      paste(colnames(x), collapse = ", ")
    )
  }

  x %>%
    mutate(
      TF_ID = as.character(.data$TF_ID),
      Gene_ID = as.character(.data$Gene_ID),
      Family = normalize_family(.data$Family),
      TF_ID_clean = normalize_id(.data$TF_ID),
      Gene_ID_clean = normalize_id(.data$Gene_ID),
      TF_ID_no_iso = str_replace(.data$TF_ID_clean, "\\.[0-9]+$", ""),
      Gene_ID_no_iso = str_replace(.data$Gene_ID_clean, "\\.[0-9]+$", "")
    ) %>%
    distinct(.data$TF_ID, .data$TF_ID_clean, .data$TF_ID_no_iso,
             .data$Gene_ID, .data$Gene_ID_clean, .data$Gene_ID_no_iso,
             .data$Family)
}

map_plant_family <- function(motif_id, plant_map) {
  motif_raw <- as.character(motif_id)
  motif_clean <- normalize_id(motif_raw)
  motif_no_iso <- str_replace(motif_clean, "\\.[0-9]+$", "")

  fam <- plant_map$Family[match(motif_raw, plant_map$TF_ID)]

  idx <- is.na(fam)
  fam[idx] <- plant_map$Family[match(motif_clean[idx], plant_map$TF_ID_clean)]

  idx <- is.na(fam)
  fam[idx] <- plant_map$Family[match(motif_no_iso[idx], plant_map$TF_ID_no_iso)]

  idx <- is.na(fam)
  fam[idx] <- plant_map$Family[match(motif_raw[idx], plant_map$Gene_ID)]

  idx <- is.na(fam)
  fam[idx] <- plant_map$Family[match(motif_clean[idx], plant_map$Gene_ID_clean)]

  idx <- is.na(fam)
  fam[idx] <- plant_map$Family[match(motif_no_iso[idx], plant_map$Gene_ID_no_iso)]

  normalize_family(fam)
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

collapse_one_family <- function(df) {
  df <- df %>%
    filter(!is.na(.data$hit_start0), !is.na(.data$hit_end0)) %>%
    arrange(.data$hit_start0, .data$hit_end0)

  if (nrow(df) == 0) return(tibble())

  starts <- df$hit_start0
  ends <- df$hit_end0

  block_id <- integer(length(starts))
  current_block <- 1L
  current_end <- ends[1]
  block_id[1] <- current_block

  if (length(starts) > 1) {
    for (i in 2:length(starts)) {
      if (starts[i] <= current_end) {
        block_id[i] <- current_block
        current_end <- max(current_end, ends[i], na.rm = TRUE)
      } else {
        current_block <- current_block + 1L
        current_end <- ends[i]
        block_id[i] <- current_block
      }
    }
  }

  df$family_block_id_local <- block_id

  df %>%
    group_by(.data$family_block_id_local) %>%
    summarise(
      window_seqnames = dplyr::first(.data$window_seqnames),
      window_start0 = dplyr::first(.data$window_start0),
      window_end0 = dplyr::first(.data$window_end0),
      side = dplyr::first(.data$side),
      geneId = dplyr::first(.data$geneId),
      tss = dplyr::first(.data$tss),
      window_strand = dplyr::first(.data$window_strand),
      block_start0 = min(.data$hit_start0, na.rm = TRUE),
      block_end0 = max(.data$hit_end0, na.rm = TRUE),
      n_raw_tfbs_hits_in_block = n(),
      motif_ids = paste(sort(unique(.data$motif_id)), collapse = ";"),
      .groups = "drop"
    )
}

make_modules_one_window <- function(df) {
  df <- df %>%
    arrange(.data$block_start0, .data$block_end0, .data$TF_family)

  if (nrow(df) < min(ngram_lengths)) return(tibble())

  out <- list()
  idx <- 1L

  for (n in ngram_lengths) {
    if (nrow(df) < n) next

    for (i in seq_len(nrow(df) - n + 1L)) {
      sub <- df[i:(i + n - 1L), , drop = FALSE]
      gaps <- sub$block_start0[-1] - sub$block_end0[-n]

      if (all(!is.na(gaps)) && all(gaps >= min_gap_bp) && all(gaps <= max_gap_bp)) {
        families <- sub$TF_family

        out[[idx]] <- tibble(
          ngram_length = n,
          module_start0 = min(sub$block_start0, na.rm = TRUE),
          module_end0 = max(sub$block_end0, na.rm = TRUE),
          gaps_bp = paste(gaps, collapse = ","),
          ordered_module = paste(families, collapse = " > "),
          unordered_multiset_module = paste(sort(families), collapse = " + "),
          unordered_presence_module = paste(sort(unique(families)), collapse = " + "),
          n_distinct_families = n_distinct(families),
          block_start0_list = paste(sub$block_start0, collapse = ","),
          block_end0_list = paste(sub$block_end0, collapse = ","),
          block_family_list = paste(families, collapse = ",")
        )
        idx <- idx + 1L
      }
    }
  }

  if (length(out) == 0) tibble() else bind_rows(out)
}

make_module_enrichment <- function(module_instances, module_col, module_type, pos_total, bg_total) {
  module_instances %>%
    transmute(
      set_name = .data$set_name,
      fasta_name = .data$fasta_name,
      ngram_length = .data$ngram_length,
      module = .data[[module_col]]
    ) %>%
    distinct(.data$set_name, .data$fasta_name, .data$ngram_length, .data$module) %>%
    count(.data$ngram_length, .data$module, .data$set_name, name = "n_windows_with_module") %>%
    pivot_wider(
      names_from = set_name,
      values_from = n_windows_with_module,
      values_fill = 0
    ) %>%
    mutate(
      positive = replace_na(.data$positive, 0L),
      background = replace_na(.data$background, 0L),
      module_type = module_type,
      positive_total = pos_total,
      background_total = bg_total,
      positive_with_module = .data$positive,
      background_with_module = .data$background,
      positive_without_module = .data$positive_total - .data$positive_with_module,
      background_without_module = .data$background_total - .data$background_with_module,
      positive_fraction = .data$positive_with_module / .data$positive_total,
      background_fraction = .data$background_with_module / .data$background_total,
      fold_positive_vs_background = (.data$positive_fraction + 1e-12) / (.data$background_fraction + 1e-12),
      fisher_p_positive_enriched = purrr::pmap_dbl(
        list(
          .data$positive_with_module,
          .data$positive_without_module,
          .data$background_with_module,
          .data$background_without_module
        ),
        ~ safe_fisher_p(..1, ..2, ..3, ..4, alternative = "greater")
      ),
      odds_ratio_positive_vs_background = purrr::pmap_dbl(
        list(
          .data$positive_with_module,
          .data$positive_without_module,
          .data$background_with_module,
          .data$background_without_module
        ),
        ~ safe_fisher_or(..1, ..2, ..3, ..4, alternative = "greater")
      ),
      passes_min_positive_support = .data$positive_with_module >= min_positive_support_for_reporting
    ) %>%
    group_by(.data$module_type, .data$ngram_length) %>%
    mutate(
      fisher_padj_BH_positive_enriched = p.adjust(.data$fisher_p_positive_enriched, method = "BH")
    ) %>%
    ungroup() %>%
    select(
      .data$module_type,
      .data$ngram_length,
      .data$module,
      .data$positive_with_module,
      .data$positive_without_module,
      .data$background_with_module,
      .data$background_without_module,
      .data$positive_total,
      .data$background_total,
      .data$positive_fraction,
      .data$background_fraction,
      .data$fold_positive_vs_background,
      .data$fisher_p_positive_enriched,
      .data$fisher_padj_BH_positive_enriched,
      .data$odds_ratio_positive_vs_background,
      .data$passes_min_positive_support
    )
}


# ============================================================
# 4) Read windows and intersected TFBS
# ============================================================

pos_windows <- read_windows_bed(pos_windows_bed)
bg_windows <- read_windows_bed(bg_windows_bed)

pos_total <- nrow(pos_windows)
bg_total <- nrow(bg_windows)

message("positive windows: ", pos_total)
message("background windows: ", bg_total)

tfbs_pos <- read_intersect(pos_intersect, "positive", tfbs_ncol)
tfbs_bg <- read_intersect(bg_intersect, "background", tfbs_ncol)

tfbs_all <- bind_rows(tfbs_pos, tfbs_bg)

readr::write_tsv(
  tfbs_all,
  file.path(table_dir, "01_intersected_known_tfbs_positive_background_raw.tsv")
)


# ============================================================
# 5) Map motif IDs to TF family
# ============================================================

plant_map <- read_plant_tf_map(family_map_file)

tfbs_family <- tfbs_all %>%
  mutate(
    motif_id_clean = normalize_id(.data$motif_id),

    # Your PlantTFDB BED name column already looks like:
    #   C2H2|MLOC_1876
    #   MYB;MYB_related|MLOC_52439
    # Therefore the TF family can be parsed directly from the part before "|".
    TF_family_from_bed = if_else(
      stringr::str_detect(.data$motif_id, "\\|"),
      stringr::str_replace(.data$motif_id, "\\|.*$", ""),
      NA_character_
    ),
    TF_family_from_bed = normalize_family(.data$TF_family_from_bed),

    TF_family_from_map = map_plant_family(.data$motif_id, plant_map),

    TF_family = dplyr::coalesce(.data$TF_family_from_bed, .data$TF_family_from_map),
    TF_family = if_else(is.na(.data$TF_family) | .data$TF_family == "", "UNKNOWN", .data$TF_family)
  )

family_map_diagnostic <- tfbs_family %>%
  count(.data$set_name, .data$TF_family, name = "n_hits") %>%
  arrange(.data$set_name, desc(.data$n_hits))

unknown_summary <- tfbs_family %>%
  summarise(
    n_hits_total = n(),
    n_hits_unknown = sum(.data$TF_family == "UNKNOWN", na.rm = TRUE),
    pct_hits_unknown = 100 * n_hits_unknown / n_hits_total,
    n_distinct_motif_ids = n_distinct(.data$motif_id),
    n_distinct_unknown_motif_ids = n_distinct(.data$motif_id[.data$TF_family == "UNKNOWN"])
  )

readr::write_tsv(
  tfbs_family,
  file.path(table_dir, "02_intersected_known_tfbs_with_tf_family.tsv")
)

readr::write_tsv(
  family_map_diagnostic,
  file.path(table_dir, "03_tf_family_raw_hit_counts.tsv")
)

readr::write_tsv(
  unknown_summary,
  file.path(table_dir, "04_tf_family_mapping_unknown_summary.tsv")
)

if (drop_unknown_family) {
  tfbs_family_used <- tfbs_family %>% filter(.data$TF_family != "UNKNOWN")
} else {
  tfbs_family_used <- tfbs_family
}


# ============================================================
# 6) Merge overlapping same-family hits into TFBS-family blocks
# ============================================================

message("Merging overlapping same-family TFBS hits into blocks...")

family_blocks <- tfbs_family_used %>%
  group_by(.data$set_name, .data$fasta_name, .data$TF_family) %>%
  group_modify(~ collapse_one_family(.x)) %>%
  ungroup() %>%
  arrange(.data$set_name, .data$fasta_name, .data$block_start0, .data$block_end0, .data$TF_family) %>%
  mutate(
    family_block_id = paste(.data$set_name, .data$fasta_name, row_number(), sep = "|")
  )

readr::write_tsv(
  family_blocks,
  file.path(table_dir, "05_tfbs_family_merged_blocks_positive_background.tsv")
)

block_summary <- family_blocks %>%
  group_by(.data$set_name) %>%
  summarise(
    n_windows_total = if_else(dplyr::first(.data$set_name) == "positive", pos_total, bg_total),
    n_windows_with_any_family_block = n_distinct(.data$fasta_name),
    n_family_blocks = n(),
    n_distinct_families = n_distinct(.data$TF_family),
    .groups = "drop"
  )

readr::write_tsv(
  block_summary,
  file.path(table_dir, "06_tfbs_family_block_summary.tsv")
)


# ============================================================
# 7) Build modules with gap_min5_max50
# ============================================================

message("Building modules with min_gap=5 bp and max_gap=50 bp...")

module_instances <- family_blocks %>%
  group_by(.data$set_name, .data$fasta_name) %>%
  group_modify(~ make_modules_one_window(.x)) %>%
  ungroup()

readr::write_tsv(
  module_instances,
  file.path(table_dir, "07_tfbs_family_module_instances_gap_min10_max50.tsv")
)

module_instance_summary <- module_instances %>%
  group_by(.data$set_name, .data$ngram_length) %>%
  summarise(
    n_module_instances = n(),
    n_windows_with_any_module = n_distinct(.data$fasta_name),
    n_ordered_modules = n_distinct(.data$ordered_module),
    n_unordered_multiset_modules = n_distinct(.data$unordered_multiset_module),
    n_unordered_presence_modules = n_distinct(.data$unordered_presence_module),
    .groups = "drop"
  )

readr::write_tsv(
  module_instance_summary,
  file.path(table_dir, "08_tfbs_family_module_instance_summary.tsv")
)


# ============================================================
# 8) Fisher positive vs all-CAGE background
# ============================================================

module_enrichment <- bind_rows(
  make_module_enrichment(module_instances, "ordered_module", "ordered", pos_total, bg_total),
  make_module_enrichment(module_instances, "unordered_multiset_module", "unordered_multiset", pos_total, bg_total),
  make_module_enrichment(module_instances, "unordered_presence_module", "unordered_presence", pos_total, bg_total)
) %>%
  arrange(.data$module_type, .data$ngram_length, .data$fisher_padj_BH_positive_enriched, desc(.data$fold_positive_vs_background))

module_enrichment_filtered <- module_enrichment %>%
  filter(.data$passes_min_positive_support) %>%
  arrange(.data$module_type, .data$ngram_length, .data$fisher_padj_BH_positive_enriched, desc(.data$fold_positive_vs_background))

readr::write_tsv(
  module_enrichment,
  file.path(table_dir, "09_module_enrichment_positive_vs_background.tsv")
)

readr::write_tsv(
  module_enrichment_filtered,
  file.path(table_dir, "10_module_enrichment_positive_vs_background_min_support.tsv")
)


# ============================================================
# 9) Family-level presence enrichment
# ============================================================

family_presence <- family_blocks %>%
  distinct(.data$set_name, .data$fasta_name, .data$TF_family)

family_enrichment <- family_presence %>%
  count(.data$TF_family, .data$set_name, name = "n_windows_with_family") %>%
  pivot_wider(
    names_from = set_name,
    values_from = n_windows_with_family,
    values_fill = 0
  ) %>%
  mutate(
    positive = replace_na(.data$positive, 0L),
    background = replace_na(.data$background, 0L),
    positive_total = pos_total,
    background_total = bg_total,
    positive_without_family = .data$positive_total - .data$positive,
    background_without_family = .data$background_total - .data$background,
    positive_fraction = .data$positive / .data$positive_total,
    background_fraction = .data$background / .data$background_total,
    fold_positive_vs_background = (.data$positive_fraction + 1e-12) / (.data$background_fraction + 1e-12),
    fisher_p_positive_enriched = purrr::pmap_dbl(
      list(.data$positive, .data$positive_without_family, .data$background, .data$background_without_family),
      ~ safe_fisher_p(..1, ..2, ..3, ..4, alternative = "greater")
    ),
    odds_ratio_positive_vs_background = purrr::pmap_dbl(
      list(.data$positive, .data$positive_without_family, .data$background, .data$background_without_family),
      ~ safe_fisher_or(..1, ..2, ..3, ..4, alternative = "greater")
    ),
    fisher_padj_BH_positive_enriched = p.adjust(.data$fisher_p_positive_enriched, method = "BH")
  ) %>%
  arrange(.data$fisher_padj_BH_positive_enriched, desc(.data$fold_positive_vs_background))

readr::write_tsv(
  family_enrichment,
  file.path(table_dir, "11_tf_family_presence_enrichment_positive_vs_background.tsv")
)


# ============================================================
# 10) Run summary
# ============================================================

run_summary <- tibble(
  min_gap_bp = min_gap_bp,
  max_gap_bp = max_gap_bp,
  ngram_lengths = paste(ngram_lengths, collapse = ","),
  positive_total_windows = pos_total,
  background_total_windows = bg_total,
  tfbs_ncol = tfbs_ncol,
  tfbs_motif_col = tfbs_motif_col,
  n_intersected_tfbs_hits = nrow(tfbs_all),
  n_tfbs_hits_with_family_used = nrow(tfbs_family_used),
  n_family_blocks = nrow(family_blocks),
  n_module_instances = nrow(module_instances),
  n_module_tests_total = nrow(module_enrichment),
  n_module_tests_min_positive_support = nrow(module_enrichment_filtered),
  min_positive_support_for_reporting = min_positive_support_for_reporting,
  drop_unknown_family = drop_unknown_family
)

readr::write_tsv(
  run_summary,
  file.path(table_dir, "00_known_tfbs_module_run_summary.tsv")
)

message("\nDONE known TFBS module enrichment")
message("Output: ", table_dir)

message("\nRun summary:")
print(run_summary)

message("\nUnknown family mapping summary:")
print(unknown_summary)

message("\nModule instance summary:")
print(module_instance_summary)

message("\nTop enriched modules with min positive support:")
print(
  module_enrichment_filtered %>%
    filter(.data$fisher_padj_BH_positive_enriched < 0.05) %>%
    group_by(.data$module_type, .data$ngram_length) %>%
    slice_head(n = 15) %>%
    ungroup()
)

message("\nTop enriched TF families:")
print(
  family_enrichment %>%
    filter(.data$fisher_padj_BH_positive_enriched < 0.05) %>%
    slice_head(n = 30)
)



# ============================================================
# Variant-specific/shared analysis of enriched known TFBS modules
# ============================================================

# ------------------------------------------------------------
# 1) Settings
# ------------------------------------------------------------

selected_module_type <- "unordered_presence"
selected_ngram_length <- 2L
selected_fdr_cutoff <- 0.05
min_any_support_for_reporting <- 5L

# Use automatic selection from the known TFBS module Fisher table.
auto_select_modules <- TRUE

# Fallback if auto_select_modules <- FALSE.
selected_modules_manual <- c(
  "ERF",
  "C2H2 + ERF",
  "ERF + MYB",
  "BBR-BPC + C2H2",
  "C2H2 + LBD",
  "ERF + LBD",
  "C2H2"
)

# ------------------------------------------------------------
# 2) Paths inherited from the known TFBS module-enrichment section
# ------------------------------------------------------------

fisher_file <- file.path(table_dir, "09_module_enrichment_positive_vs_background.tsv")
module_instance_file <- file.path(table_dir, "07_tfbs_family_module_instances_gap_min10_max50.tsv")

if (!file.exists(fisher_file)) {
  stop("Cannot find module enrichment table: ", fisher_file)
}

if (!file.exists(module_instance_file)) {
  stop("Cannot find module instance table: ", module_instance_file)
}

out_dir <- table_dir

message("Using module enrichment table: ", fisher_file)
message("Using module instance table: ", module_instance_file)


# ------------------------------------------------------------
# 3) Helper functions
# ------------------------------------------------------------

read_fasta_headers <- function(path) {
  hdr <- readLines(path, warn = FALSE)
  hdr <- hdr[startsWith(hdr, ">")]
  hdr <- sub("^>", "", hdr)

  tibble(sequence_name_header = hdr) %>%
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
      tss = suppressWarnings(as.integer(str_match(.data$fasta_name, "TSS(\\d+)")[, 2]))
    )
}

parse_positive_fasta_name <- function(x) {
  tibble(fasta_name = x) %>%
    mutate(
      geneId = str_match(.data$fasta_name, "^[^|]+\\|([^|]+)\\|")[, 2],
      side = str_match(.data$fasta_name, "^[^|]+\\|[^|]+\\|([^|]+)\\|")[, 2],
      tss = suppressWarnings(as.integer(str_match(.data$fasta_name, "TSS(\\d+)")[, 2]))
    )
}

safe_binom_p <- function(variant_specific, any_total) {
  if (is.na(variant_specific) || is.na(any_total) || any_total <= 0) return(NA_real_)
  tryCatch(
    stats::binom.test(
      x = variant_specific,
      n = any_total,
      p = 0.5,
      alternative = "greater"
    )$p.value,
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

# ------------------------------------------------------------
# 4) Gene universe: all genes with >=1 long and >=1 short
# ------------------------------------------------------------

pos_headers <- read_fasta_headers(pos_fasta)

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

write_tsv(gene_side_window_counts, file.path(out_dir, "12_gene_long_short_window_counts.tsv"))
write_tsv(gene_set_summary, file.path(out_dir, "13_gene_set_summary.tsv"))

# ------------------------------------------------------------
# 5) Select the enriched modules to test
# ------------------------------------------------------------

fisher_tab <- read_tsv(fisher_file, show_col_types = FALSE)

if (auto_select_modules) {
  selected_modules <- fisher_tab %>%
    filter(
      .data$module_type == selected_module_type,
      .data$ngram_length == selected_ngram_length,
      .data$fisher_padj_BH_positive_enriched < selected_fdr_cutoff
    )

  if ("passes_min_positive_support" %in% colnames(selected_modules)) {
    selected_modules <- selected_modules %>% filter(.data$passes_min_positive_support)
  }

  selected_modules <- selected_modules %>%
    arrange(.data$fisher_padj_BH_positive_enriched, desc(.data$fold_positive_vs_background))
} else {
  selected_modules <- fisher_tab %>%
    filter(
      .data$module_type == selected_module_type,
      .data$ngram_length == selected_ngram_length,
      .data$module %in% selected_modules_manual
    ) %>%
    arrange(.data$fisher_padj_BH_positive_enriched, desc(.data$fold_positive_vs_background))
}

if (nrow(selected_modules) == 0) stop("No selected modules found.")

write_tsv(selected_modules, file.path(out_dir, "14_selected_modules_for_variant_specific_test.tsv"))

# ------------------------------------------------------------
# 6) Read module instances and annotate positive windows
# ------------------------------------------------------------

module_instances <- read_tsv(module_instance_file, show_col_types = FALSE)

module_instances_annot <- module_instances %>%
  filter(.data$set_name == "positive") %>%
  left_join(
    parse_positive_fasta_name(unique(module_instances$fasta_name)),
    by = "fasta_name"
  )

module_instances_long <- module_instances_annot %>%
  select(
    .data$set_name,
    .data$fasta_name,
    .data$geneId,
    .data$side,
    .data$tss,
    .data$ngram_length,
    .data$ordered_module,
    .data$unordered_multiset_module,
    .data$unordered_presence_module
  ) %>%
  pivot_longer(
    cols = c("ordered_module", "unordered_multiset_module", "unordered_presence_module"),
    names_to = "module_type_raw",
    values_to = "module"
  ) %>%
  mutate(
    module_type = case_when(
      .data$module_type_raw == "ordered_module" ~ "ordered",
      .data$module_type_raw == "unordered_multiset_module" ~ "unordered_multiset",
      .data$module_type_raw == "unordered_presence_module" ~ "unordered_presence",
      TRUE ~ .data$module_type_raw
    )
  ) %>%
  inner_join(
    selected_modules %>% select(.data$module_type, .data$ngram_length, .data$module),
    by = c("module_type", "ngram_length", "module")
  ) %>%
  distinct(
    .data$module_type,
    .data$ngram_length,
    .data$module,
    .data$geneId,
    .data$side,
    .data$fasta_name,
    .data$tss
  )

write_tsv(module_instances_long, file.path(out_dir, "15_selected_module_instances_positive_windows_annotated.tsv"))

# ------------------------------------------------------------
# 7) Gene-level module presence: any long vs any short
# ------------------------------------------------------------

module_gene_side_presence <- module_instances_long %>%
  filter(!is.na(.data$geneId), .data$side %in% c("long", "short")) %>%
  count(
    .data$module_type,
    .data$ngram_length,
    .data$module,
    .data$geneId,
    .data$side,
    name = "n_windows_with_module_on_side"
  ) %>%
  pivot_wider(
    names_from = side,
    values_from = n_windows_with_module_on_side,
    values_fill = 0,
    names_prefix = "n_module_windows_"
  ) %>%
  mutate(
    n_module_windows_long = replace_na(.data$n_module_windows_long, 0L),
    n_module_windows_short = replace_na(.data$n_module_windows_short, 0L),
    long = .data$n_module_windows_long > 0,
    short = .data$n_module_windows_short > 0
  )

module_gene_grid <- tidyr::crossing(
  selected_modules %>% select(.data$module_type, .data$ngram_length, .data$module),
  long_short_gene_set %>% select(.data$geneId, .data$n_long, .data$n_short, .data$gene_pair_class)
)

gene_level_detail <- module_gene_grid %>%
  left_join(
    module_gene_side_presence,
    by = c("module_type", "ngram_length", "module", "geneId")
  ) %>%
  mutate(
    n_module_windows_long = replace_na(.data$n_module_windows_long, 0L),
    n_module_windows_short = replace_na(.data$n_module_windows_short, 0L),
    long = replace_na(.data$long, FALSE),
    short = replace_na(.data$short, FALSE),
    module_presence_class = case_when(
      .data$long & .data$short ~ "shared_both_long_and_short",
      .data$long & !.data$short ~ "long_only",
      !.data$long & .data$short ~ "short_only",
      TRUE ~ "neither"
    ),
    variant_specific = xor(.data$long, .data$short),
    shared = .data$long & .data$short,
    any_variant = .data$long | .data$short
  ) %>%
  arrange(.data$module_type, .data$ngram_length, .data$module, .data$geneId)

write_tsv(gene_level_detail, file.path(out_dir, "16_module_variant_specific_shared_by_gene_detail.tsv"))

# ------------------------------------------------------------
# 8) Summary statistics
# ------------------------------------------------------------

variant_summary <- gene_level_detail %>%
  group_by(.data$module_type, .data$ngram_length, .data$module) %>%
  summarise(
    gene_total_with_at_least_one_long_and_short = n_distinct(.data$geneId),

    shared_both_long_and_short = sum(.data$shared, na.rm = TRUE),
    long_only = sum(.data$long & !.data$short, na.rm = TRUE),
    short_only = sum(!.data$long & .data$short, na.rm = TRUE),
    variant_specific = sum(.data$variant_specific, na.rm = TRUE),
    neither = sum(!.data$long & !.data$short, na.rm = TRUE),

    genes_with_module_in_long = sum(.data$long, na.rm = TRUE),
    genes_with_module_in_short = sum(.data$short, na.rm = TRUE),
    genes_with_module_in_any_variant = sum(.data$any_variant, na.rm = TRUE),

    fraction_shared_among_genes_with_module_any_variant =
      if_else(.data$genes_with_module_in_any_variant > 0,
              .data$shared_both_long_and_short / .data$genes_with_module_in_any_variant,
              NA_real_),

    fraction_variant_specific_among_genes_with_module_any_variant =
      if_else(.data$genes_with_module_in_any_variant > 0,
              .data$variant_specific / .data$genes_with_module_in_any_variant,
              NA_real_),

    expected_shared_under_independent_marginals =
      hyper_expected_shared(
        .data$genes_with_module_in_long,
        .data$genes_with_module_in_short,
        .data$gene_total_with_at_least_one_long_and_short
      ),

    sd_shared_under_independent_marginals =
      hyper_sd_shared(
        .data$genes_with_module_in_long,
        .data$genes_with_module_in_short,
        .data$gene_total_with_at_least_one_long_and_short
      ),

    shared_observed_minus_expected =
      .data$shared_both_long_and_short - .data$expected_shared_under_independent_marginals,

    shared_observed_over_expected =
      if_else(.data$expected_shared_under_independent_marginals > 0,
              .data$shared_both_long_and_short / .data$expected_shared_under_independent_marginals,
              NA_real_),

    p_shared_enriched_hypergeom =
      hyper_p_shared_enriched(
        .data$shared_both_long_and_short,
        .data$genes_with_module_in_long,
        .data$genes_with_module_in_short,
        .data$gene_total_with_at_least_one_long_and_short
      ),

    p_variant_specific_enriched_hypergeom =
      hyper_p_variant_specific_enriched(
        .data$shared_both_long_and_short,
        .data$genes_with_module_in_long,
        .data$genes_with_module_in_short,
        .data$gene_total_with_at_least_one_long_and_short
      ),

    binom_p_variant_specific_gt_shared =
      safe_binom_p(
        .data$variant_specific,
        .data$genes_with_module_in_any_variant
      ),

    mcnemar_p_long_vs_short_asymmetry =
      safe_mcnemar_p(.data$long_only, .data$short_only),

    passes_min_any_support =
      .data$genes_with_module_in_any_variant >= min_any_support_for_reporting,

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
  left_join(
    selected_modules %>%
      select(
        .data$module_type,
        .data$ngram_length,
        .data$module,
        .data$positive_with_module,
        .data$background_with_module,
        .data$fold_positive_vs_background,
        .data$fisher_padj_BH_positive_enriched
      ),
    by = c("module_type", "ngram_length", "module")
  ) %>%
  arrange(.data$fisher_padj_BH_positive_enriched, desc(.data$fold_positive_vs_background))

write_tsv(variant_summary, file.path(out_dir, "17_module_variant_specific_shared_summary.tsv"))

reporting_table <- variant_summary %>%
  filter(.data$passes_min_any_support) %>%
  select(
    .data$module_type,
    .data$ngram_length,
    .data$module,

    .data$positive_with_module,
    .data$background_with_module,
    .data$fold_positive_vs_background,
    .data$fisher_padj_BH_positive_enriched,

    .data$gene_total_with_at_least_one_long_and_short,
    .data$genes_with_module_in_any_variant,
    .data$shared_both_long_and_short,
    .data$variant_specific,
    .data$long_only,
    .data$short_only,
    .data$neither,

    .data$fraction_shared_among_genes_with_module_any_variant,
    .data$fraction_variant_specific_among_genes_with_module_any_variant,

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
  ) %>%
  arrange(.data$fisher_padj_BH_positive_enriched, desc(.data$fold_positive_vs_background))

write_tsv(reporting_table, file.path(out_dir, "18_module_variant_specific_shared_reporting_table.tsv"))

run_summary <- tibble(
  input_known_tfbs_table_dir = table_dir,
  output_dir = out_dir,
  selected_module_type = selected_module_type,
  selected_ngram_length = selected_ngram_length,
  selected_fdr_cutoff = selected_fdr_cutoff,
  auto_select_modules = auto_select_modules,
  n_selected_modules = nrow(selected_modules),
  n_positive_windows_total = nrow(pos_headers),
  n_unique_genes_in_positive_headers = n_distinct(pos_headers$geneId[!is.na(pos_headers$geneId)]),
  n_genes_with_at_least_one_long_and_one_short_window = nrow(long_short_gene_set),
  min_any_support_for_reporting = min_any_support_for_reporting
)

write_tsv(run_summary, file.path(out_dir, "19_module_variant_specific_run_summary.tsv"))

message("\nDONE module variant-specific/shared analysis")
message("Output: ", out_dir)

message("\nRun summary:")
print(run_summary)

message("\nSelected modules:")
print(
  selected_modules %>%
    select(
      module_type,
      ngram_length,
      module,
      positive_with_module,
      background_with_module,
      fold_positive_vs_background,
      fisher_padj_BH_positive_enriched
    )
)

message("\nVariant-specific/shared reporting table:")
print(reporting_table)
