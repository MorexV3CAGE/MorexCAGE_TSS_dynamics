#looking for barley shifts
#---------------------------------------------------------------------------
# use packages
library(GenomicFeatures)
library(BSgenome)
library(rtracklayer)
library(ggplot2)
library(ChIPseeker)
library(dplyr)
library(BSgenome.MorexV3.Gatersleben)

#First load datasets from CAGEr that you want to analyze
df_8DAP <- read.delim("/path/to/CAGE/promoters_8DAP.bed", header = T, sep = "\t", quote = "\"", fill = TRUE)
df_24DAP <- read.delim("/path/to/CAGE/promoters_24DAP.bed", header = T, sep = "\t", quote = "\"", fill = TRUE)
df_4DAG <- read.delim("/path/to/CAGE/promoters_4DAG.bed", header = T, sep = "\t", quote = "\"", fill = TRUE)
df_INF <- read.delim("/path/to/CAGE/promoters_INF.bed", header = T, sep = "\t", quote = "\"", fill = TRUE)
#columns should have: seqnames, start, end, geneID, score, strand, dominantTSS, annotation, tpm.dominant_ctss

#load gff as a txdb
txdb <-makeTxDbFromGFF("/path/to/annotation/Hv_Morex.pgsb.Jul2020.gff3", format="gff3")
txdb

outdir_exon_boundry <- "/path/to/outdir/exon_intron_boundry/"

annotated_8DAP <- process_and_annotate(df_8DAP)
annotated_24DAP <- process_and_annotate(df_24DAP)
annotated_4DAG <- process_and_annotate(df_4DAG)
annotated_INF <- process_and_annotate(df_INF)

prepared_8DAP <- prepare_data(annotated_8DAP, outdir_exon_boundry, "8DAP")
prepared_24DAP <- prepare_data(annotated_24DAP, outdir_exon_boundry, "24DAP")
prepared_4DAG <- prepare_data(annotated_4DAG, outdir_exon_boundry, "4DAG")
prepared_INF <- prepare_data(annotated_INF, outdir_exon_boundry, "INF")

# Prepare lists for shifts
df_list_8DAP <- list(
  "24DAP" = prepared_24DAP,
  "4DAG" = prepared_4DAG,
  "INF" = prepared_INF
)

df_list_24DAP <- list(
  "8DAP" = prepared_8DAP,
  "4DAG" = prepared_4DAG,
  "INF" = prepared_INF
)

df_list_4DAG <- list(
  "8DAP" = prepared_8DAP,
  "24DAP" = prepared_24DAP,
  "INF" = prepared_INF
)

df_list_INF <- list(
  "8DAP" = prepared_8DAP,
  "24DAP" = prepared_24DAP,
  "4DAG" = prepared_4DAG
)

# Apply processing to each sample
results_list <- list(
  "sample_8DAP" = process_sample_shifts(prepared_8DAP, df_list_8DAP, "8DAP"),
  "sample_24DAP" = process_sample_shifts(prepared_24DAP, df_list_24DAP, "24DAP"),
  "sample_4DAG" = process_sample_shifts(prepared_4DAG, df_list_4DAG, "4DAG"),
  "sample_INF" = process_sample_shifts(prepared_INF, df_list_INF, "INF")
)

shifts_8DAP_24DAP <- results_list_new$sample_8DAP$shifts_with_other_samples$shift_with_24DAP %>% filter(shift_bp <= 5000)

shifts_8DAP_4DAG <- results_list_new$sample_8DAP$shifts_with_other_samples$shift_with_4DAG %>% filter(shift_bp <= 5000)

shifts_24DAP_4DAG <- results_list_new$sample_24DAP$shifts_with_other_samples$shift_with_4DAG %>% filter(shift_bp <= 5000)

shifts_INF_8DAP <- results_list_new$sample_INF$shifts_with_other_samples$shift_with_8DAP %>% filter(shift_bp <= 5000)

shifts_INF_24DAP <- results_list_new$sample_INF$shifts_with_other_samples$shift_with_24DAP %>% filter(shift_bp <= 5000)

shifts_INF_4DAG <- results_list_new$sample_INF$shifts_with_other_samples$shift_with_4DAG %>% filter(shift_bp <= 5000)

# Categorize shifts for each comparison
shifts_8DAP_24DAP_categorized <- categorize_shifts(shifts_8DAP_24DAP)
shifts_8DAP_4DAG_categorized <- categorize_shifts(shifts_8DAP_4DAG)
shifts_24DAP_4DAG_categorized <- categorize_shifts(shifts_24DAP_4DAG)
shifts_INF_8DAP_categorized <- categorize_shifts(shifts_INF_8DAP)
shifts_INF_24DAP_categorized <- categorize_shifts(shifts_INF_24DAP)
shifts_INF_4DAG_categorized <- categorize_shifts(shifts_INF_4DAG)

# Define output directory
output_dir <- "/path/to/outdir/"

# Export all categorized shifts with signal validation
export_categorized_shifts(shifts_8DAP_24DAP_categorized, "8DAPx24DAP", output_dir, "8DAP", "24DAP")

export_categorized_shifts(shifts_8DAP_4DAG_categorized, "8DAPx4DAG", output_dir, "8DAP", "4DAG")

export_categorized_shifts(shifts_24DAP_4DAG_categorized, "24DAPx4DAG", output_dir, "24DAP", "4DAG")

export_categorized_shifts(shifts_INF_8DAP_categorized, "INF_8DAP", output_dir, "INF", "8DAP")

export_categorized_shifts(shifts_INF_24DAP_categorized, "INF_24DAP", output_dir, "INF", "24DAP")

export_categorized_shifts(shifts_INF_4DAG_categorized, "INF_4DAG", output_dir, "INF", "24DAP")


#---------------------------------------------------------------------------
#USED FUNCTIONS
#---------------------------------------------------------------------------
# Convert data frame to GRanges
df_to_granges <- function(df) {
  gr <- GRanges(
    seqnames = df$seqnames,
    ranges = IRanges(start = df$start, end = df$end), 
    strand = df$strand,
    mcols = df[, !colnames(df) %in% c("seqnames", "start", "end", "strand")] # Retain metadata columns
  )
  return(gr)
}

# Annotate peaks and process the results
annotate_and_filter <- function(gr, original_columns) {
  # Annotate using Chipseeker
  annotated <- annotatePeak(gr, tssRegion = c(-500, 100), TxDb = txdb, overlap = "all")

  annotated_df <- as.data.frame(annotated)
  
  # Remove the 'mcols.' prefix from column names
  colnames(annotated_df) <- sub("^mcols\\.", "", colnames(annotated_df))
  
  # Check to remove any wrong annotation
  annotated_df <- annotated_df[annotated_df$geneID == annotated_df$geneId,]
  
  # Map geneStrand (1/2) to strand symbols (+/-)
  annotated_df$geneStrand <- ifelse(annotated_df$geneStrand == 1, "+", "-")
  
  # Remove rows with antisense transcripts
  filtered_df <- annotated_df[annotated_df$strand == annotated_df$geneStrand, ]
  
  # Retain only the columns from the original data frame
  final_df <- filtered_df[, original_columns, drop = FALSE]
  
  return(final_df)
}

# Annotate the data and remove antisense transcripts
process_and_annotate <- function(df) {
  # Store the original column names
  original_columns <- colnames(df)
  gr <- df_to_granges(df)
  annotated_df <- annotate_and_filter(gr, original_columns)
  return(annotated_df)
}

# Function to remove exon/intron boundries from the analysis
check_tss_near_exon_starts <- function(shift_df, txdb, output_dir, sample_name) {
  # Extract exons and group by gene to identify the first exon
  exons_by_gene <- exonsBy(txdb, by = "gene")
  
  remove_first_exon <- function(exons_by_gene) {
    all_exons <- unlist(exons_by_gene)
    gene_ids <- names(all_exons)
    strands <- strand(all_exons)
    
    # Remove 1st exon since we are interested in exon beginnings from 2nd exon
    is_first_exon <- !duplicated(gene_ids)
    is_last_exon <- !duplicated(gene_ids, fromLast = TRUE)
    to_remove <- (is_first_exon & strands == "+") | (is_last_exon & strands == "-")
    return(all_exons[!to_remove])
  }
  
  # Apply the function
  non_first_exons <- remove_first_exon(exons_by_gene)
  exons_df <- as.data.frame(non_first_exons, row.names = c(1:length(non_first_exons)))
  
  # Adjust exon start positions based on strand, for the remaining exons
  exons_df <- exons_df %>%
    mutate(
      start = if_else(strand == "+", start, end)  # Use start for + strand, end for - strand
    )
  
  # Convert to GRanges
  exon_starts <- GRanges(
    seqnames = exons_df$seqnames,
    ranges = IRanges(start = exons_df$start, end = exons_df$start),
    strand = exons_df$strand
  )
  
  # Expand exon starts by ±5bp (5bp before and after start)
  expanded_exon_starts <- exon_starts %>%
    resize(width = 11, fix = "center")
  
  # Create GRanges for dominantTSS from the shift dataframe
  tss_ranges <- GRanges(
    seqnames = shift_df$seqnames,
    ranges = IRanges(start = shift_df$dominantTSS, end = shift_df$dominantTSS),
    strand = shift_df$strand
  )
  
  # Check overlaps between dominantTSS and expanded exon starts
  overlaps <- findOverlaps(tss_ranges, expanded_exon_starts)
  
  # Add a new column to indicate proximity to exon starts (excluding first exon)
  shift_df <- shift_df %>%
    mutate(
      near_exon_start = row_number() %in% queryHits(overlaps)
    )
  
  # Create the output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Save the exon/intron shifts to a file
  exon_intron_shifts <- shift_df %>% filter(near_exon_start == TRUE)
  write.table(exon_intron_shifts, paste0(output_dir, "exon_intron_shifts", sample_name, ".bed"), 
              sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  return(shift_df)
}

# Function to prepare the data for the shift analysis
prepare_data <- function(combined_df, outdir_exon_boundry, sample_name){
  # Assuming you have a TxDb object (`txdb`)
  combined_df <- check_tss_near_exon_starts(combined_df, txdb, outdir_exon_boundry, sample_name)
  combined_df <- combined_df[combined_df$near_exon_start == FALSE,]
  
  #select dominantTSS
  combined_df <- combined_df %>%
    group_by(geneID) %>%
    slice_max(order_by = abs(tpm.dominant_ctss), n = 1) %>%
    ungroup()
  return(combined_df)
}

# Updated calculate_shifts function
calculate_shifts <- function(secondary_df, primary_df) {
  # Determine the correct column names for dominantTSS in each dataframe
  dominantTSS_secondary <- if ("dominantTSS" %in% colnames(secondary_df)) "dominantTSS" else "dominant"
  dominantTSS_primary <- if ("dominantTSS" %in% colnames(primary_df)) "dominantTSS" else "dominant"
  
  # Perform inner join on the geneID column
  combined_df <- inner_join(secondary_df, primary_df, by = "geneID", suffix = c("_secondary", "_primary"))
  
  # Calculate the shift in base pairs
  combined_df <- combined_df %>%
    mutate(shift_bp = abs(combined_df$dominantTSS_secondary - combined_df$dominantTSS_primary))
  
  return(combined_df)
}


# Update process_sample_shifts
process_sample_shifts <- function(combined_df, other_primary_dfs, sample_name) {
  # Step 1: Create shift dataframes
  shift_dfs <- lapply(names(other_primary_dfs), function(other_sample_name) {
    # Calculate shifts between selected secondary and primary data from other samples
    shifts <- calculate_shifts(combined_df, other_primary_dfs[[other_sample_name]])
    return(shifts)
  })
  
  # Step 2: Combine results in a list
  results <- list(
    filtered_secondary_genes = combined_df,
    shifts_with_other_samples = shift_dfs
  )
  
  # Name the shifts with other samples
  names(results$shifts_with_other_samples) <- paste("shift_with_", names(other_primary_dfs), sep = "")
  
  return(results)
}

# Function to categorize shifts
categorize_shifts <- function(shift_df) {
  # Filter out shifts where shift_bp == 0
  shift_df <- shift_df %>% filter(shift_bp != 0)
  
  # Categorize shifts into "true" and "short"
  true_shifts <- shift_df %>% filter(abs(shift_bp) >= 100)
  short_shifts <- shift_df %>% filter(abs(shift_bp) < 100)
  
  # Return a list of categorized shifts
  return(list(
    true_shifts = true_shifts,
    short_shifts = short_shifts
  ))
}

# Function to export categorized shifts
export_categorized_shifts <- function(categorized_shifts, prefix, output_dir, sample1, sample2) {
  if (nrow(categorized_shifts$true_shifts) > 0) {
    true_shifts <- categorized_shifts$true_shifts
    
    export_shift_dfs_to_bed_single(
      true_shifts,
      file.path(output_dir, paste0(prefix, "_true_shifts.bed")), 
      sample1, sample2)
  }
  
  if (nrow(categorized_shifts$short_shifts) > 0) {
    short_shifts <- categorized_shifts$short_shifts
    
    export_shift_dfs_to_bed_single(
      short_shifts,
      file.path(output_dir, paste0(prefix, "_short_shifts.bed")), 
      sample1, sample2)
  }
}

# Define the export function
export_shift_dfs_to_bed_single <- function(shift_df, output_file, sample1, sample2) {
  # Check and adjust start and end values
  shift_df <- shift_df %>%
    mutate(
      start = pmin(dominantTSS_secondary, dominantTSS_primary),
      end = pmax(dominantTSS_secondary, dominantTSS_primary),
      fromto = ifelse(strand_secondary == "+", 
                      ifelse(start == dominantTSS_secondary, 
                             paste0("close_", sample2, "_far_", sample1), 
                             paste0("close_", sample1, "_far_", sample2)), 
                      ifelse(start == dominantTSS_secondary, 
                             paste0("close_", sample1, "_far_", sample2), 
                             paste0("close_", sample2, "_far_", sample1)))
    )
  
  # Select relevant columns for BED format (chromosome, start, end, geneID, shift_bp)
  bed_df <- shift_df %>% 
    select(chrom = seqnames_secondary, start = start, end = end, name = geneID, score = shift_bp, strand = strand_secondary, annotation1 = annotation_secondary, annotation2 = annotation_primary, TPM1 = tpm.dominant_ctss_secondary, TPM2 = tpm.dominant_ctss_primary, logFC = log2fc, fromto = fromto)
  
  # Load gene bed file
  filter_genes <- read.table("/path/to/reference/genes.bed", header = FALSE, 
                             col.names = c("chr", "start", "end", "name", "score", "strand"))
  
  # Convert to GRanges objects for overlap detection
  bed_gr <- GRanges(seqnames = bed_df$chrom,
                    ranges = IRanges(start = bed_df$start, end = bed_df$end))
  
  filter_gr <- GRanges(seqnames = filter_genes$chr,
                       ranges = IRanges(start = filter_genes$start, end = filter_genes$end),
                       strand = filter_genes$strand)
  
  # Find overlaps between bed_gr and filter_gr
  overlaps <- findOverlaps(bed_gr, filter_gr)
  
  # Create a data frame of overlaps
  overlap_df <- as.data.frame(overlaps)
  
  # Count overlaps per row in bed_df (i.e., how many times each row from bed_df overlaps)
  overlap_count <- table(overlap_df$queryHits)
  
  # Filter out rows from bed_df that overlap with at least 2 filter_bed_data rows
  filtered_bed_df <- bed_df[!rownames(bed_df) %in% names(overlap_count[overlap_count >= 2]), ]
  
  # Remove geneIDs that appear more than once in filtered_bed_df
  gene_counts <- table(filtered_bed_df$name)
  unique_filtered_bed_df <- filtered_bed_df[!filtered_bed_df$name %in% names(gene_counts[gene_counts > 1]), ]
  
  # Write to BED file
  write.table(unique_filtered_bed_df, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  message("Export complete. Files are saved in ", output_file)
}
