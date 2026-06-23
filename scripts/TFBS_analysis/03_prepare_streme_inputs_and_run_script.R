# ============================================================
# 03) Prepare STREME input sequences and run script
#
# Goal:
#   Prepare positive/background promoter-window BED files for de novo
#   motif discovery, and generate a shell script that extracts FASTA
#   sequences and runs STREME from MEME Suite.
#
# Biological design:
#   positive set:
#     shifted promoter variants
#
#   background set:
#     all CAGE-defined promoter/TSS entries, optionally excluding exact
#     positive matches
#
#   region:
#     + strand: [TSS - 150, TSS - 50)
#     - strand: [TSS + 50, TSS + 150)
#
# Run:
#   Rscript 03_prepare_streme_inputs_and_run_script.R
#   bash results/streme_motif_discovery/commands/run_streme.sh
# ============================================================



# ============================================================
# 0) Libraries
# ============================================================

library(dplyr)
library(readr)
library(stringr)
library(tidyr)

if (requireNamespace("conflicted", quietly = TRUE)) {
  conflicted::conflicts_prefer(base::setdiff)
  conflicted::conflicts_prefer(base::intersect)
  conflicted::conflicts_prefer(base::union)

  conflicted::conflicts_prefer(dplyr::filter)
  conflicted::conflicts_prefer(dplyr::select)
  conflicted::conflicts_prefer(dplyr::mutate)
  conflicted::conflicts_prefer(dplyr::arrange)
  conflicted::conflicts_prefer(dplyr::count)
  conflicted::conflicts_prefer(dplyr::distinct)
  conflicted::conflicts_prefer(dplyr::summarise)
  conflicted::conflicts_prefer(dplyr::anti_join)
  conflicted::conflicts_prefer(dplyr::left_join)
}


# ============================================================
# 1) Settings
# ============================================================

project_dir <- Sys.getenv("PROJECT_DIR", unset = ".")
results_dir <- Sys.getenv("RESULTS_DIR", unset = file.path(project_dir, "results"))

# Input tables.
cage_file <- Sys.getenv(
  "CAGE_PROMOTER_TABLE",
  unset = file.path(project_dir, "input", "cage", "all_cage_promoters.bed")
)

shift_file <- Sys.getenv(
  "SHIFT_PROMOTER_TABLE",
  unset = file.path(project_dir, "input", "shift_promoters", "shift_promoter_variants.bed")
)

# Output directory for STREME input preparation and STREME run.
out_dir <- Sys.getenv(
  "STREME_DISCOVERY_OUTPUT_DIR",
  unset = file.path(results_dir, "streme_motif_discovery")
)

table_dir <- file.path(out_dir, "tables")
bed_dir   <- file.path(out_dir, "bed")
fasta_dir <- file.path(out_dir, "fasta")
cmd_dir   <- file.path(out_dir, "commands")
log_dir   <- file.path(out_dir, "logs")

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(bed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fasta_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cmd_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

# Window definition.
upstream_start_bp <- 50L
upstream_end_bp   <- 150L

# Keep only long/short rows if a side column exists.
# Set FALSE if the shift table contains only the desired variants.
keep_only_long_short <- TRUE

# Background design.
exclude_positives_from_background <- TRUE

# If TRUE, background uses only CAGE rows annotated as Promoter.
# Default FALSE = all CAGE consensus entries.
background_annotation_promoter_only <- FALSE

# Paths used in the generated shell script.
# They are resolved when this R script is run; edit config_example.sh or set
# these environment variables before running if needed.
genome_fasta <- Sys.getenv(
  "GENOME_FASTA",
  unset = "/path/to/reference_genome.fa"
)

meme_sif <- Sys.getenv(
  "MEME_SIF",
  unset = "/path/to/meme.sif"
)

bedtools_cmd <- Sys.getenv(
  "BEDTOOLS_CMD",
  unset = "bedtools"
)

# STREME parameters in the generated shell script.
streme_minw <- as.integer(Sys.getenv("STREME_MINW", unset = "6"))
streme_maxw <- as.integer(Sys.getenv("STREME_MAXW", unset = "15"))
streme_seed <- as.integer(Sys.getenv("STREME_SEED", unset = "123"))
streme_threshold <- Sys.getenv("STREME_THRESH", unset = "0.05")
streme_use_evalue <- as.logical(Sys.getenv("STREME_USE_EVALUE", unset = "TRUE"))


# ============================================================
# 2) Helper functions
# ============================================================

find_first_existing_col <- function(df, candidates, label, required = TRUE) {
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) == 0) {
    if (required) {
      stop(
        "Cannot find column for ", label, ". Tried: ",
        paste(candidates, collapse = ", ")
      )
    } else {
      return(NA_character_)
    }
  }
  hit[1]
}

clean_tss <- function(x) {
  x <- as.character(x)
  x <- stringr::str_remove(x, "^TSSpos_")
  as.integer(x)
}

standardize_strand <- function(x) {
  x <- as.character(x)
  dplyr::case_when(
    x %in% c("+", "plus", "1", "forward", "F") ~ "+",
    x %in% c("-", "minus", "-1", "reverse", "R") ~ "-",
    TRUE ~ NA_character_
  )
}

make_unique_name <- function(x) {
  ave(
    x,
    x,
    FUN = function(z) {
      if (length(z) == 1) {
        z
      } else {
        paste0(z, "__dup", seq_along(z))
      }
    }
  )
}

make_bed_name_positive <- function(df) {
  base <- paste(
    "POS",
    df$geneId,
    df$side,
    paste0("TSS", df$tss_pos),
    sep = "|"
  )

  if ("shift_name" %in% colnames(df)) {
    base <- paste(base, df$shift_name, sep = "|")
  }

  base <- gsub("[[:space:]]+", "_", base)
  base <- gsub("[;]", "_", base)
  make_unique_name(base)
}

make_bed_name_background <- function(df) {
  base <- paste(
    "BG",
    df$geneId,
    paste0("TSS", df$dominantTSS),
    df$type,
    df$annotation,
    sep = "|"
  )

  base <- gsub("[[:space:]]+", "_", base)
  base <- gsub("[;]", "_", base)
  make_unique_name(base)
}

make_upstream_bed_from_shift <- function(df) {
  df %>%
    dplyr::mutate(
      bed_start = dplyr::case_when(
        .data$strand == "+" ~ .data$tss_pos - upstream_end_bp,
        .data$strand == "-" ~ .data$tss_pos + upstream_start_bp,
        TRUE ~ NA_integer_
      ),
      bed_end = dplyr::case_when(
        .data$strand == "+" ~ .data$tss_pos - upstream_start_bp,
        .data$strand == "-" ~ .data$tss_pos + upstream_end_bp,
        TRUE ~ NA_integer_
      ),
      bed_start = pmax(0L, as.integer(.data$bed_start)),
      bed_end = as.integer(.data$bed_end),
      width = .data$bed_end - .data$bed_start,
      name = make_bed_name_positive(.),
      score = 0L
    ) %>%
    dplyr::filter(
      !is.na(.data$seqnames),
      !is.na(.data$bed_start),
      !is.na(.data$bed_end),
      .data$bed_end > .data$bed_start,
      .data$strand %in% c("+", "-"),
      .data$width == upstream_end_bp - upstream_start_bp
    ) %>%
    dplyr::transmute(
      seqnames = as.character(.data$seqnames),
      start = as.integer(.data$bed_start),
      end = as.integer(.data$bed_end),
      name = as.character(.data$name),
      score = as.integer(.data$score),
      strand = as.character(.data$strand)
    )
}

make_upstream_bed_from_cage <- function(df) {
  df %>%
    dplyr::mutate(
      bed_start = dplyr::case_when(
        .data$strand == "+" ~ .data$dominantTSS - upstream_end_bp,
        .data$strand == "-" ~ .data$dominantTSS + upstream_start_bp,
        TRUE ~ NA_integer_
      ),
      bed_end = dplyr::case_when(
        .data$strand == "+" ~ .data$dominantTSS - upstream_start_bp,
        .data$strand == "-" ~ .data$dominantTSS + upstream_end_bp,
        TRUE ~ NA_integer_
      ),
      bed_start = pmax(0L, as.integer(.data$bed_start)),
      bed_end = as.integer(.data$bed_end),
      width = .data$bed_end - .data$bed_start,
      name = make_bed_name_background(.),
      score = 0L
    ) %>%
    dplyr::filter(
      !is.na(.data$seqnames),
      !is.na(.data$bed_start),
      !is.na(.data$bed_end),
      .data$bed_end > .data$bed_start,
      .data$strand %in% c("+", "-"),
      .data$width == upstream_end_bp - upstream_start_bp
    ) %>%
    dplyr::transmute(
      seqnames = as.character(.data$seqnames),
      start = as.integer(.data$bed_start),
      end = as.integer(.data$bed_end),
      name = as.character(.data$name),
      score = as.integer(.data$score),
      strand = as.character(.data$strand)
    )
}

write_bed6 <- function(bed, path) {
  readr::write_tsv(
    bed,
    path,
    col_names = FALSE
  )
}


# ============================================================
# 3) Load inputs
# ============================================================

if (!file.exists(cage_file)) {
  stop("Missing CAGE file: ", cage_file)
}

if (!file.exists(shift_file)) {
  stop("Missing shift file: ", shift_file)
}

cage <- readr::read_tsv(cage_file, show_col_types = FALSE)
shift <- readr::read_tsv(shift_file, show_col_types = FALSE)

message("Loaded CAGE rows: ", nrow(cage))
message("Loaded shift rows: ", nrow(shift))

required_cage <- c("seqnames", "geneId", "strand", "dominantTSS", "annotation", "type")
missing_cage <- base::setdiff(required_cage, colnames(cage))
if (length(missing_cage) > 0) {
  stop("CAGE file is missing columns: ", paste(missing_cage, collapse = ", "))
}

shift_chr_col <- find_first_existing_col(
  shift,
  c("seqnames", "chr", "chrom", "chromosome"),
  "shift chromosome"
)

shift_pos_col <- find_first_existing_col(
  shift,
  c("pos", "TSS_pos", "TSS_position", "dominantTSS", "tss", "TSS"),
  "shift TSS position"
)

shift_strand_col <- find_first_existing_col(
  shift,
  c("gene_strand", "strand", "TSS_strand"),
  "shift strand",
  required = FALSE
)

if (!"geneId" %in% colnames(shift)) {
  stop("Shift file is missing required column: geneId")
}

if (!"side" %in% colnames(shift)) {
  shift$side <- "shift"
}

if (!"shift_name" %in% colnames(shift)) {
  shift$shift_name <- paste0("shift_", seq_len(nrow(shift)))
}

shift_std <- shift %>%
  dplyr::mutate(
    seqnames = as.character(.data[[shift_chr_col]]),
    geneId = as.character(.data$geneId),
    tss_pos = clean_tss(.data[[shift_pos_col]]),
    side = as.character(.data$side),
    shift_name = as.character(.data$shift_name)
  )

if (keep_only_long_short && "side" %in% colnames(shift_std)) {
  shift_std <- shift_std %>%
    dplyr::filter(.data$side %in% c("long", "short"))
}

cage_std <- cage %>%
  dplyr::mutate(
    seqnames = as.character(.data$seqnames),
    geneId = as.character(.data$geneId),
    strand = standardize_strand(.data$strand),
    dominantTSS = as.integer(.data$dominantTSS),
    annotation = as.character(.data$annotation),
    type = as.character(.data$type)
  ) %>%
  dplyr::filter(.data$strand %in% c("+", "-"))

# ============================================================
# 4) Add strand to shift variants
# ============================================================

if (!is.na(shift_strand_col)) {
  shift_std <- shift_std %>%
    dplyr::mutate(
      strand = standardize_strand(shift[[shift_strand_col]]),
      strand_source = dplyr::if_else(!is.na(.data$strand), "shift_file", NA_character_)
    )
} else {
  shift_std <- shift_std %>%
    dplyr::mutate(
      strand = NA_character_,
      strand_source = NA_character_
    )
}

# Try exact CAGE geneId + TSS match for missing strand.
cage_strand_exact <- cage_std %>%
  dplyr::select(
    .data$geneId,
    tss_pos = .data$dominantTSS,
    cage_strand = .data$strand
  ) %>%
  dplyr::distinct()

shift_std <- shift_std %>%
  dplyr::left_join(
    cage_strand_exact,
    by = c("geneId", "tss_pos")
  ) %>%
  dplyr::mutate(
    strand = dplyr::if_else(
      is.na(.data$strand) & !is.na(.data$cage_strand),
      .data$cage_strand,
      .data$strand
    ),
    strand_source = dplyr::if_else(
      is.na(.data$strand_source) & !is.na(.data$cage_strand),
      "CAGE_exact_geneId_TSS",
      .data$strand_source
    )
  ) %>%
  dplyr::select(-.data$cage_strand)

# Infer strand from long/short order within gene for remaining missing strand.
strand_inferred_by_gene <- shift_std %>%
  dplyr::filter(.data$side %in% c("long", "short")) %>%
  dplyr::group_by(.data$geneId) %>%
  dplyr::summarise(
    median_long_tss = median(.data$tss_pos[.data$side == "long"], na.rm = TRUE),
    median_short_tss = median(.data$tss_pos[.data$side == "short"], na.rm = TRUE),
    n_long = sum(.data$side == "long", na.rm = TRUE),
    n_short = sum(.data$side == "short", na.rm = TRUE),
    inferred_strand = dplyr::case_when(
      n_long > 0 & n_short > 0 & median_long_tss < median_short_tss ~ "+",
      n_long > 0 & n_short > 0 & median_long_tss > median_short_tss ~ "-",
      TRUE ~ NA_character_
    ),
    .groups = "drop"
  ) %>%
  dplyr::select(.data$geneId, .data$inferred_strand)

shift_std <- shift_std %>%
  dplyr::left_join(strand_inferred_by_gene, by = "geneId") %>%
  dplyr::mutate(
    strand = dplyr::if_else(
      is.na(.data$strand) & !is.na(.data$inferred_strand),
      .data$inferred_strand,
      .data$strand
    ),
    strand_source = dplyr::if_else(
      is.na(.data$strand_source) & !is.na(.data$inferred_strand),
      "inferred_long_short_order_by_gene",
      .data$strand_source
    )
  ) %>%
  dplyr::select(-.data$inferred_strand)

# Remove exact duplicated shift rows only, not biological variants.
shift_std <- shift_std %>%
  dplyr::distinct(
    .data$shift_name,
    .data$geneId,
    .data$side,
    .data$seqnames,
    .data$tss_pos,
    .keep_all = TRUE
  )

strand_diag <- shift_std %>%
  dplyr::count(.data$strand_source, .data$strand, name = "n_variants") %>%
  dplyr::arrange(.data$strand_source, .data$strand)

readr::write_tsv(
  strand_diag,
  file.path(table_dir, "01_positive_strand_source_diagnostic.tsv")
)

unresolved_strand <- shift_std %>%
  dplyr::filter(is.na(.data$strand) | !.data$strand %in% c("+", "-"))

if (nrow(unresolved_strand) > 0) {
  readr::write_tsv(
    unresolved_strand,
    file.path(table_dir, "01B_positive_unresolved_strand_rows.tsv")
  )
  warning(
    "Some positive shift variants have unresolved strand and will be dropped: ",
    nrow(unresolved_strand),
    ". See 01B_positive_unresolved_strand_rows.tsv"
  )
}

positive_shift <- shift_std %>%
  dplyr::filter(.data$strand %in% c("+", "-"))

readr::write_tsv(
  positive_shift,
  file.path(table_dir, "02_positive_shift_promoter_variants.tsv")
)


# ============================================================
# 5) Define background from all CAGE consensus rows
# ============================================================

background_cage <- cage_std

if (background_annotation_promoter_only) {
  background_cage <- background_cage %>%
    dplyr::filter(.data$annotation == "Promoter")
}

if (exclude_positives_from_background) {
  positive_keys <- positive_shift %>%
    dplyr::transmute(
      bg_key = paste(.data$seqnames, .data$geneId, .data$tss_pos, .data$strand, sep = "|")
    ) %>%
    dplyr::distinct()

  background_cage <- background_cage %>%
    dplyr::mutate(
      bg_key = paste(.data$seqnames, .data$geneId, .data$dominantTSS, .data$strand, sep = "|")
    ) %>%
    dplyr::anti_join(positive_keys, by = "bg_key") %>%
    dplyr::select(-.data$bg_key)
}

readr::write_tsv(
  background_cage,
  file.path(table_dir, "03_background_cage_promoter_rows.tsv")
)


# ============================================================
# 6) Create BED files
# ============================================================

positive_bed <- make_upstream_bed_from_shift(positive_shift)
background_bed <- make_upstream_bed_from_cage(background_cage)

positive_bed_file <- file.path(bed_dir, "positive_shift_promoter_windows_50_150bp_up.bed")
background_bed_file <- file.path(bed_dir, "background_all_cage_promoter_windows_50_150bp_up.bed")

write_bed6(positive_bed, positive_bed_file)
write_bed6(background_bed, background_bed_file)

message("\nPositive shift rows after side filter/dedup: ", nrow(shift_std))
message("Positive shift rows with resolved strand: ", nrow(positive_shift))
message("Positive BED rows: ", nrow(positive_bed))
message("Background CAGE rows: ", nrow(background_cage))
message("Background BED rows: ", nrow(background_bed))

summary_df <- tibble::tibble(
  cage_file = cage_file,
  shift_file = shift_file,
  shift_chr_column_used = shift_chr_col,
  shift_tss_column_used = shift_pos_col,
  shift_strand_column_used = ifelse(is.na(shift_strand_col), "", shift_strand_col),
  shift_rows_loaded = nrow(shift),
  positive_shift_rows_after_side_filter_and_dedup = nrow(shift_std),
  positive_shift_rows_with_resolved_strand = nrow(positive_shift),
  positive_bed_rows = nrow(positive_bed),
  background_cage_rows = nrow(background_cage),
  background_bed_rows = nrow(background_bed),
  upstream_start_bp = upstream_start_bp,
  upstream_end_bp = upstream_end_bp,
  keep_only_long_short = keep_only_long_short,
  exclude_positives_from_background = exclude_positives_from_background,
  background_annotation_promoter_only = background_annotation_promoter_only,
  positive_bed_file = positive_bed_file,
  background_bed_file = background_bed_file
)

readr::write_tsv(
  summary_df,
  file.path(table_dir, "00_streme_input_preparation_summary.tsv")
)

if (nrow(positive_bed) < 2) {
  stop("Positive BED has fewer than 2 sequences.")
}

if (nrow(background_bed) < 2) {
  stop("Background BED has fewer than 2 sequences.")
}


# ============================================================
# 7) Write shell script for FASTA extraction and STREME
# ============================================================

positive_fasta <- file.path(fasta_dir, "positive_shift_promoter_windows_50_150bp_up.fa")
background_fasta <- file.path(fasta_dir, "background_all_cage_promoter_windows_50_150bp_up.fa")
streme <- file.path(out_dir, "streme")

run_script <- file.path(cmd_dir, "run_streme.sh")

shell_lines <- c(
  "#!/usr/bin/env bash",
  "set -euo pipefail",
  "",
  "# ============================================================",
  "# STREME motif discovery",
  "# ============================================================",
  "",
  "# These values were written by 03_prepare_streme_inputs_and_run_script.R.",
  "# They can still be overridden by environment variables before running.",
  paste0("GENOME_FASTA=\"", genome_fasta, "\""),
  paste0("MEME_SIF=\"", meme_sif, "\""),
  paste0("BEDTOOLS_CMD=\"", bedtools_cmd, "\""),
  "",
  paste0("POS_BED=\"", positive_bed_file, "\""),
  paste0("BG_BED=\"", background_bed_file, "\""),
  paste0("POS_FASTA=\"", positive_fasta, "\""),
  paste0("BG_FASTA=\"", background_fasta, "\""),
  paste0("STREME_OUT=\"", streme_out, "\""),
  paste0("LOG_DIR=\"", log_dir, "\""),
  "",
  "mkdir -p \"$(dirname \"$POS_FASTA\")\" \"$(dirname \"$BG_FASTA\")\" \"$STREME_OUT\" \"$LOG_DIR\"",
  "",
  "if [[ ! -f \"$GENOME_FASTA\" ]]; then",
  "  echo \"ERROR: missing GENOME_FASTA: $GENOME_FASTA\"",
  "  exit 1",
  "fi",
  "",
  "if [[ ! -f \"$MEME_SIF\" ]]; then",
  "  echo \"ERROR: missing MEME_SIF: $MEME_SIF\"",
  "  exit 1",
  "fi",
  "",
  "echo '[1/4] Extract positive FASTA'",
  "\"$BEDTOOLS_CMD\" getfasta \",
  "  -fi \"$GENOME_FASTA\" \",
  "  -bed \"$POS_BED\" \",
  "  -s \",
  "  -name \",
  "  -fo \"$POS_FASTA\" \",
  "  > \"$LOG_DIR/bedtools_positive.log\" 2>&1",
  "",
  "echo '[2/4] Extract background FASTA'",
  "\"$BEDTOOLS_CMD\" getfasta \",
  "  -fi \"$GENOME_FASTA\" \",
  "  -bed \"$BG_BED\" \",
  "  -s \",
  "  -name \",
  "  -fo \"$BG_FASTA\" \",
  "  > \"$LOG_DIR/bedtools_background.log\" 2>&1",
  "",
  "echo '[3/4] Count FASTA sequences'",
  "echo -n 'positive sequences: '",
  "grep -c '^>' \"$POS_FASTA\"",
  "echo -n 'background sequences: '",
  "grep -c '^>' \"$BG_FASTA\"",
  "",
  "echo '[4/4] Run STREME via Singularity MEME Suite'",
  "singularity exec \"$MEME_SIF\" streme \",
  "  --dna \",
  "  --p \"$POS_FASTA\" \",
  "  --n \"$BG_FASTA\" \",
  "  --oc \"$STREME_OUT\" \",
  paste0("  --minw ", streme_minw, " \"),
  paste0("  --maxw ", streme_maxw, " \"),
  paste0("  --seed ", streme_seed, " \"),
  "  --align right \",
  if (isTRUE(streme_use_evalue)) "  --evalue \" else NULL,
  paste0("  --thresh ", streme_threshold, " \"),
  "  > \"$LOG_DIR/streme.log\" 2>&1",
  "",
  "echo 'DONE'",
  "echo \"STREME HTML: $STREME_OUT/streme.html\"",
  "echo \"STREME motifs: $STREME_OUT/streme.txt\"",
  "if [[ -f \"$STREME_OUT/streme.txt\" && ! -f \"$STREME_OUT/streme_FINAL.txt\" ]]; then",
  "  cp \"$STREME_OUT/streme.txt\" \"$STREME_OUT/streme_FINAL.txt\"",
  "fi"
)

writeLines(shell_lines, con = run_script)
Sys.chmod(run_script, mode = "0755")

message("\nGenerated shell script:")
message(run_script)
message("\nEdit GENOME_FASTA and MEME_SIF in the shell script, then run:")
message("bash ", run_script)

message("\nDONE STREME input preparation")
print(summary_df)
