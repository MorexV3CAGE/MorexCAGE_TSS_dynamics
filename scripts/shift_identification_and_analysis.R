## ============================================================================
## Barley CAGE TSS-shift pipeline
## SECTION A: shift identification
## SECTION B: differential usage of detected shifts with DESeq2
## ============================================================================

## ---- CONFIG -----------------------------------------------------------------

input_dir   <- "data/promoters"
input_files <- c(
  "8DAP"  = file.path(input_dir, "promoters_8DAP.bed"),
  "24DAP" = file.path(input_dir, "promoters_24DAP.bed"),
  "4DAG"  = file.path(input_dir, "promoters_4DAG.bed"),
  "INF"   = file.path(input_dir, "promoters_INF.bed"))
gff_file    <- "annotation/Hv_Morex.pgsb.Jul2020.gff3"
gene_bed    <- "reference/genes.bed"
outdir_exon <- "results/exon_intron_boundary/"
output_dir  <- "results/shifts/"
max_shift_bp <- 5000

## DESeq2
deseq_dir  <- "results/deseq2"
reference_cond <- "D8"

suppressPackageStartupMessages({
  library(GenomicFeatures)   # TxDb, exonsBy, findOverlaps (GenomicRanges)
  library(ChIPseeker)        # annotatePeak
  library(dplyr)
  library(DESeq2)
})

bed_raw <- lapply(input_files,
                  read.delim, header = TRUE, sep = "\t", quote = "\"", fill = TRUE)
## expected columns: seqnames, start, end, geneID, score, strand, dominantTSS, annotation, tpm.dominant_ctss

txdb <- makeTxDbFromGFF(gff_file, format = "gff3")

## ============================================================================
## SHIFT IDENTIFICATION
## ============================================================================

## ---- Filter to remove shifts between multiple close genes ------
genes_filter <- read.table(gene_bed, header = FALSE,
  col.names = c("chr", "start", "end", "name", "score", "strand"))
genes_filter_gr <- GRanges(
  seqnames = genes_filter$chr,
  ranges   = IRanges(start = genes_filter$start, end = genes_filter$end),
  strand   = genes_filter$strand)

## ---- Filter for exon/intron boundary signals --------
non_first_exon_starts <- local({
  all_exons <- unlist(exonsBy(txdb, by = "gene"))
  gids      <- names(all_exons)
  strds     <- strand(all_exons)
  drop <- (!duplicated(gids)                  & strds == "+") |
          (!duplicated(gids, fromLast = TRUE) & strds == "-")
  keep <- all_exons[!drop]
  kdf  <- as.data.frame(keep)
  kdf$start <- ifelse(kdf$strand == "+", kdf$start, kdf$end)
  resize(GRanges(seqnames = kdf$seqnames,
                 IRanges(start = kdf$start, end = kdf$start),
                 strand = kdf$strand),
         width = 11, fix = "center")
})

## ---- Annotate and filter antisense (if not already filtered) -----------
annotate_and_filter <- function(df) {
  gr <- GRanges(seqnames = df$seqnames,
                ranges   = IRanges(start = df$start, end = df$end),
                strand   = df$strand,
                mcols    = df[, setdiff(colnames(df),
                               c("seqnames", "start", "end", "strand"))])
  a  <- annotatePeak(gr, tssRegion = c(-500, 100), TxDb = txdb,
                     overlap = "all")
  ad <- as.data.frame(a)
  colnames(ad) <- sub("^mcols\\.", "", colnames(ad))
  ad <- ad[ad$geneID == ad$geneId, ]
  ad$geneStrand <- ifelse(ad$geneStrand == 1, "+", "-")
  ## drop antisense, restore original columns
  ad[ad$strand == ad$geneStrand, colnames(df), drop = FALSE]
}

## ---- Prepare samples: exon/intron boundary filter + dominantTSS pick ----------
prepare_sample <- function(df, sample_name) {
  tss_ranges <- GRanges(seqnames = df$seqnames,
                        ranges   = IRanges(start = df$dominantTSS,
                                           end   = df$dominantTSS),
                        strand   = df$strand)
  df$near_exon_start <- seq_len(nrow(df)) %in%
    queryHits(findOverlaps(tss_ranges, non_first_exon_starts))

  write.table(dplyr::filter(df, near_exon_start),
              file.path(outdir_exon,
                        paste0("exon_intron_shifts", sample_name, ".bed")),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

  df <- df[!df$near_exon_start, ]
  df %>%
    group_by(geneID) %>%
    slice_max(order_by = abs(tpm.dominant_ctss), n = 1) %>%
    ungroup()
}

annotated <- lapply(bed_raw, annotate_and_filter)
prepared  <- Map(prepare_sample, annotated, names(annotated))

## ---- Shift detection between sample pairs ----------------------------------
calculate_shifts <- function(secondary_df, primary_df) {
  inner_join(secondary_df, primary_df,
             by = "geneID", suffix = c("_secondary", "_primary")) %>%
    mutate(shift_bp = abs(dominantTSS_secondary - dominantTSS_primary))
}

sample_names <- names(prepared)
results_list <- setNames(
  lapply(sample_names, function(sn) {
    list(filtered_secondary_genes = prepared[[sn]],
         shifts_with_other_samples = setNames(
           lapply(sample_names[sample_names != sn],
                  function(other) calculate_shifts(prepared[[sn]],
                                                   prepared[[other]])),
           paste0("shift_with_", sample_names[sample_names != sn])))
  }),
  paste0("sample_", sample_names))

comparisons <- list(
  c(d1 = "8DAP",  d2 = "24DAP", label = "8DAPx24DAP"),
  c(d1 = "8DAP",  d2 = "4DAG",  label = "8DAPx4DAG"),
  c(d1 = "24DAP", d2 = "4DAG",  label = "24DAPx4DAG"),
  c(d1 = "INF",   d2 = "8DAP",  label = "INFx8DAP"),
  c(d1 = "INF",   d2 = "24DAP", label = "INFx24DAP"),
  c(d1 = "INF",   d2 = "4DAG",  label = "INFx4DAG"))

## Build the six shift tables, filtered for length (shift_bp <= 5000)
shift_tables <- lapply(comparisons, function(cmp) {
  results_list[[paste0("sample_", cmp[["d1"]])]]$shifts_with_other_samples[[
    paste0("shift_with_", cmp[["d2"]])]] %>%
    filter(shift_bp <= max_shift_bp, shift_bp != 0)
})
names(shift_tables) <- vapply(comparisons, `[[`, "", "label")

categorize_shifts <- function(shift_df) {
  list(true_shifts  = shift_df %>% filter(abs(shift_bp) >= 100),
       short_shifts = shift_df %>% filter(abs(shift_bp) < 100))
}

## ---- Export -----------------------------------------------------
export_shifts <- function(shift_df, output_file, sample1, sample2) {
  if (nrow(shift_df) == 0) return(invisible(NULL))

  shift_df <- shift_df %>%
    mutate(start = pmin(dominantTSS_secondary, dominantTSS_primary),
           end   = pmax(dominantTSS_secondary, dominantTSS_primary),
           fromto = ifelse(strand_secondary == "+",
                    ifelse(start == dominantTSS_secondary,
                           paste0("close_", sample2, "_far_", sample1),
                           paste0("close_", sample1, "_far_", sample2)),
                    ifelse(start == dominantTSS_secondary,
                           paste0("close_", sample1, "_far_", sample2),
                           paste0("close_", sample2, "_far_", sample1))))

  bed_df <- shift_df %>%
    select(chrom = seqnames_secondary, start, end, name = geneID,
           score = shift_bp, strand = strand_secondary,
           annotation1 = annotation_secondary,
           annotation2 = annotation_primary,
           TPM1 = tpm.dominant_ctss_secondary,
           TPM2 = tpm.dominant_ctss_primary,
           fromto)

  ## Filter the close genes here
  bed_gr <- GRanges(seqnames = bed_df$chrom,
                    ranges   = IRanges(start = bed_df$start, end = bed_df$end))
  ov_tab <- as.data.frame(findOverlaps(bed_gr, genes_filter_gr))
  ov_cnt <- table(ov_tab$queryHits)
  bed_df <- bed_df[!rownames(bed_df) %in% names(ov_cnt[ov_cnt >= 2]), ]

  gene_cnt <- table(bed_df$name)
  bed_df   <- bed_df[!bed_df$name %in% names(gene_cnt[gene_cnt > 1]), ]

  write.table(bed_df, output_file, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE)
  message("Export complete: ", output_file)
}

for (cmp in comparisons) {
  cat <- categorize_shifts(shift_tables[[cmp[["label"]]]])
  export_shifts(cat$true_shifts,
                file.path(output_dir,
                          paste0(cmp[["label"]], "_true_shifts.bed")),
                cmp[["d1"]], cmp[["d2"]])
  export_shifts(cat$short_shifts,
                file.path(output_dir,
                          paste0(cmp[["label"]], "_short_shifts.bed")),
                cmp[["d1"]], cmp[["d2"]])
}

## ============================================================================
## DIFFERENTIAL USAGE OF DETECTED SHIFTS (DESeq2)
## ============================================================================
## INPUTS
##   1) *_true_shifts.bed files written by export_categorized_shifts()
##        (one per condition comparison, e.g. 24DAPx4DAG_true_shifts.bed)
##   2) tc from CAGEr analysis - consensusClustersGR(object)
##   3) tagcount_coords : per-CTSS tag counts - CTSStagCountTable(object)

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(GenomicRanges)
})

tagcount_file   <- "data/ctss_tag_counts.bed"

## ---- 1. Load the detected shifts ---------------------------------------------
shift_files <- list.files(
  output_dir, pattern = "true", full.names = TRUE)

shift_list <- lapply(shift_files, read.delim, header = FALSE) %>%
  setNames(
    ## turn "8DAPx24DAP_true_shifts.bed" into "8DAP_24DAP"
    sub("x",   "_",  basename(shift_files)) %>%
      sub("_true_shifts.bed", "", .)
  )

cond_map <- c(
  "8DAP"  = "D8",
  "24DAP" = "D24",
  "4DAG"  = "D4",
  "INF"   = "INF")

## ---- 2. Sum CTSS tag counts over promoters (tc) ------------------------------
tagcount_coords <- read.delim(tagcount_file, header = TRUE)
## columns expected: seqnames, pos, strand, <one column per sample>

tag_gr <- GRanges(
  seqnames = tagcount_coords$seqnames,
  ranges   = IRanges(start = tagcount_coords$pos, end = tagcount_coords$pos),
  strand   = tagcount_coords$strand
)

mcols(tag_gr) <- tagcount_coords[, setdiff(names(tagcount_coords),
                                           c("seqnames", "pos", "strand"))]

## The promoter table = the consensus clusters from CAGEr.
promoters_gr <- tc

## Which CTSS positions fall inside which promoter interval.
ov <- findOverlaps(promoters_gr, tag_gr)

## For every sample column: sum the CTSS counts per promoter
sample_cols <- setdiff(names(mcols(tag_gr)), "strand")   # e.g. D4.A, D4.B, ...
for (s in sample_cols) {
  sums <- tapply(mcols(tag_gr)[[s]][subjectHits(ov)],
                 queryHits(ov),
                 sum)
  mcols(promoters_gr)[[paste0(s, "_sum")]] <- 0L
  mcols(promoters_gr)[[paste0(s, "_sum")]][as.integer(names(sums))] <- sums
}

## Reduce the GRanges to a plain data frame: coordinates + all *_sum columns.
promoters_sum_df <- as.data.frame(promoters_gr) %>%
  select(chr = seqnames, start, end, ends_with("_sum"))

## ---- 3. Attach promoter tag counts to every detected shift ------------------
## For each shift row, the "close" endpoint is column V2 and the "far" endpoint
## is column V3 (start/end positions in the bed file). For every *_sum column
## we add two columns: "<col>_sum_V2" (counts of the promoter interval that
## contains V2) and "<col>_sum_V3" (same for V3). This tells us how strongly
## each of the two promoters of the shift is expressed in every replicate.
process_shift_df <- function(shift_df, promoters_sum_df) {
  sum_cols <- grep("_sum$", colnames(promoters_sum_df), value = TRUE)

  ## Row-wise point-in-interval lookup, one column at a time (kept from the
  ## original script: simple and correct, just slow for very large tables).
  for (colname in sum_cols) {
    shift_df[[paste0(colname, "_V2")]] <- apply(shift_df, 1, function(row) {
      vals <- promoters_sum_df[[colname]][
        promoters_sum_df$chr   == row["V1"] &
        promoters_sum_df$start <= as.numeric(row["V2"]) &
        promoters_sum_df$end   >= as.numeric(row["V2"])]
      if (length(vals) == 0) return(0L)
      sum(vals)
    })
    shift_df[[paste0(colname, "_V3")]] <- apply(shift_df, 1, function(row) {
      vals <- promoters_sum_df[[colname]][
        promoters_sum_df$chr   == row["V1"] &
        promoters_sum_df$start <= as.numeric(row["V3"]) &
        promoters_sum_df$end   >= as.numeric(row["V3"])]
      if (length(vals) == 0) return(0L)
      sum(vals)
    })
  }

  ## ---- tag_close / tag_far --------------------------------------------------
  strand_is_plus <- shift_df$V6 == "+"

  ## Parse "close_<cond>_far_<cond>" once for the whole table (vectorised).
  tag_id <- str_match(shift_df$V12, "^close_(.+?)_far_(.+)$")
  stopifnot("unexpected fromto tag in column V12" = !any(is.na(tag_id)))
  close_cond <- tag_id[, 2]
  far_cond   <- tag_id[, 3]
  stopifnot("unknown condition in fromto tag - check cond_map" =
      all(close_cond %in% names(cond_map)) && all(far_cond %in% names(cond_map)))

  tag_close <- tag_far <- numeric(nrow(shift_df))
  for (sc in sum_cols) {
    col_prefix  <- sub("\\..*$", "", sc)     # e.g. "D4.A_sum" -> "D4"
    cond_of_col <- names(cond_map)[cond_map == col_prefix]  # "D4" -> "4DAG"
    if (length(cond_of_col) == 0) next
    is_close_col <- cond_of_col == close_cond
    is_far_col   <- cond_of_col == far_cond
    if (!any(is_close_col | is_far_col)) next
    use_v3 <- (is_close_col &  strand_is_plus) |
              (is_far_col   & !strand_is_plus)
    val <- ifelse(use_v3,
                  shift_df[[paste0(sc, "_V3")]],
                  shift_df[[paste0(sc, "_V2")]])
    tag_close <- tag_close + val * is_close_col
    tag_far   <- tag_far   + val * is_far_col
  }
  shift_df$tag_close <- tag_close
  shift_df$tag_far   <- tag_far

  shift_df
}

## Apply to each comparison in the list.
shift_list_with_tags_TC <- lapply(shift_list, process_shift_df,
                                     promoters_sum_df)
saveRDS(shift_list_with_tags_TC, shift_rds)

## ---- 4. Build the DESeq2 count matrix ----------------------------------------
build_matrix <- function(input_data) {
  assign_vals <- function(df, i, smp, rep, role, strand) {
    coln <- sprintf("%s.%d_sum", smp, rep)
    flip <- (role == "close" && strand == "-") ||
            (role == "far"   && strand == "+")
    sfx  <- if (flip) c("_V2_V3", "_V3") else c("_V3", "_V2_V3")
    v1 <- suppressWarnings(as.numeric(df[[paste0(coln, sfx[1])]][i]))
    v2 <- suppressWarnings(as.numeric(df[[paste0(coln, sfx[2])]][i]))
    val <- max(c(v1, v2), na.rm = TRUE)
    if (is.infinite(val) || is.na(val)) 0 else val
  }

  rows <- list(); all_cols <- character(0)
  for (nm in names(input_data)) {
    df <- input_data[[nm]]
    df <- df[grepl("close_", df$V12) & grepl("_far_", df$V12), , drop = FALSE]
    if (nrow(df) == 0) next
    meta <- sub("close_", "", df$V12)
    for (i in seq_len(nrow(df))) {
      smp_c <- sub("_far_.*$", "", meta[i])
      smp_f <- sub(".*_far_",  "", meta[i])
      strand <- df$V6[i]
      vals <- c()
      for (role in c("close", "far")) {
        smp <- if (role == "close") smp_c else smp_f
        for (rep in c(1, 2)) {
          coln <- sprintf("%s.%d_sum", smp, rep)
          vals[coln] <- assign_vals(df, i, smp, rep, role, strand)
        }
      }
      all_cols <- union(all_cols, names(vals))
      rname <- paste(df$V4[i], df$V1[i], df$V2[i], df$V3[i], sep = "_")
      rows[[length(rows) + 1L]] <- data.frame(row = rname,
                                              t(as.data.frame(vals)))
    }
  }
  stopifnot("no usable rows: check close_/far_ tags" = length(rows) > 0)
  mat <- do.call(rbind, rows)
  for (cc in setdiff(all_cols, names(mat)[-1])) mat[[cc]] <- 0
  agg <- aggregate(. ~ row, data = mat, FUN = mean)
  m <- as.matrix(agg[, sort(names(agg)[-1]), drop = FALSE])
  rownames(m) <- agg$row
  m
}

## input
input_data <- if (exists("shift_list_with_tags_TC")) shift_list_with_tags_TC
              else readRDS(file.path(output_dir, "shift_list_with_tags_TC.rds"))

count_mat <- build_matrix(input_data)
write.csv(count_mat, file.path(deseq_dir, "shift_pair_counts.csv"))
message("Matrix: ", nrow(count_mat), " shift pairs x ", ncol(count_mat), " samples")

sample_condition <- sub("\\..*$", "", colnames(count_mat))
stopifnot("reference condition not found among samples" =
          reference_cond %in% sample_condition)

coldata <- data.frame(
  row.names = colnames(count_mat),
  condition = relevel(factor(sample_condition), ref = reference_cond))

dds <- DESeqDataSetFromMatrix(countData = round(count_mat),
                              colData   = coldata,
                              design    = ~ condition)
## Optional pre-filter: dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)

conds <- levels(coldata$condition)
for (cmp in combn(conds, 2, simplify = FALSE)) {
  name <- paste(rev(cmp), collapse = "_vs_")   # e.g. "INF_vs_D8"
  res  <- results(dds, contrast = c("condition", cmp[2], cmp[1]))
  res  <- as.data.frame(res[order(res$pvalue), ])
  write.table(res,
    file   = file.path(deseq_dir, paste0("DESeq2_", name, ".txt")),
    sep = "\t", quote = FALSE, row.names = TRUE)
}
message("Done. Differential results -> ", normalizePath(deseq_dir))