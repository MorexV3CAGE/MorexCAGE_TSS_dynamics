## ============================================================================
## CAGE analysis with CAGEr
## Based on the official CAGEr vignette:
## https://bioconductor.org/packages/release/bioc/vignettes/CAGEr/inst/doc/CAGEexp.html
## ============================================================================

## ---- CONFIG -----------------------------------------------------------------

genome_name <- "BSgenome.MorexV3"   # BSgenome package with your genome
bam_dir     <- "data/bam"           # folder containing input BAM files
gff_file    <- "annotation/genome_annotation.gff3"
out_dir     <- "results"

samples <- data.frame(
  bam   = c("groupA_rep1.bam", "groupA_rep2.bam",
            "groupB_rep1.bam", "groupB_rep2.bam"),
  group = c("groupA", "groupA", "groupB", "groupB")
)

## getCTSS()
seq_quality_threshold <- 10
map_quality_threshold <- 20
remove_first_g       <- TRUE
correct_systematic_g <- FALSE

## normaliseTagCount()
norm_fit_range <- c(10, 10000)
norm_alpha     <- 1.05   # power-law exponent: check plotReverseCumulatives()!
norm_T         <- 10^6

## clusterCTSS()
cluster_threshold    <- 0.1
cluster_nr_pass      <- 1
cluster_max_dist     <- 100
cluster_keep_singles <- 5

## quantilePositions() / aggregateTagClusters()
q_low            <- 0.1
q_up             <- 0.9
tc_tpm_threshold <- 1
tc_max_dist      <- 100

## Parallelisation (MulticoreParam is not available on Windows)
use_multicore <- TRUE
nr_cores      <- 4

## ---- Packages ---------------------------------------------------------------

suppressPackageStartupMessages({
  library(CAGEr)
  library(GenomicFeatures)
  library(ChIPseeker)
  library(DESeq2)
})

if (!genome_name %in% .packages()) library(genome_name, character.only = TRUE)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

if (use_multicore && requireNamespace("BiocParallel", quietly = TRUE)) {
  BiocParallel::register(BiocParallel::MulticoreParam(workers = nr_cores))
}

## Helper: save a ggplot-returning function to PDF
.save_plot <- function(expr, filename, width = 7, height = 5) {
  pdf(file.path(out_dir, filename), width = width, height = height)
  p <- expr
  if (!is.null(p)) print(p)
  dev.off()
}

## ---- Load data and build the CTSS object ------------------------------------

input_files   <- file.path(bam_dir, samples$bam)
sample_labels <- sub("\\.bam$", "", basename(samples$bam))

ce <- CAGEexp(genomeName     = genome_name,
              inputFiles     = input_files,
              inputFilesType = "bam",
              sampleLabels   = sample_labels)

ce <- getCTSS(ce,
              sequencingQualityThreshold = seq_quality_threshold,
              mappingQualityThreshold    = map_quality_threshold,
              removeFirstG               = remove_first_g,
              correctSystematicG         = correct_systematic_g,
              useMulticore               = use_multicore,
              nrCores                    = nr_cores)

librarySizes(ce)
CTSStagCountSE(ce)   # peek at the count matrix

## ---- Merge replicates ---------------------------------------------------------
## Optional
group_levels <- unique(samples$group)
ce_merged <- mergeSamples(ce,
                          mergeIndex         = match(samples$group, group_levels),
                          mergedSampleLabels = group_levels)


## ---- Normalization ---------------------------------------------------
## Choose fitInRange / alpha so the fitted power law matches your data.

.save_plot(plotReverseCumulatives(ce, fitInRange = norm_fit_range),
           "reverse_cumulatives.pdf")

ce <- normalizeTagCount(ce, method = "powerLaw",
                        fitInRange = norm_fit_range,
                        alpha = norm_alpha, T = norm_T)

## ---- CTSS clustering and interquantile widths --------------------------------

ce <- filterLowExpCTSS(ce, thresholdIsTpm = TRUE, nrPassThreshold = cluster_nr_pass, threshold = cluster_threshold)

ce <- distclu(ce, maxDist = cluster_max_dist, keepSingletonsAbove = cluster_keep_singles)

ce <- cumulativeCTSSdistribution(ce, clusters = "tagClusters",
                                 useMulticore = use_multicore)
ce <- quantilePositions(ce, clusters = "tagClusters",
                        qLow = q_low, qUp = q_up)

## saveRDS(ce, file.path(out_dir, "ce_clustered.rds"))


## ---- Aggregate tag clusters across samples -----------------------------------

ce <- aggregateTagClusters(ce,
                           tpmThreshold = tc_tpm_threshold,
                           qLow = q_low, qUp = q_up,
                           maxDist = tc_max_dist)

tc <- consensusClustersGR(ce)

## ---- Annotation and rDNA filter ---------------------------------------------------

tc_anno <- annotatePeak(tc, tssRegion=c(-500, 100),
                         TxDb=txdb, overlap = "all")

tc_df <- as.data.frame(tc_anno)

## masking regions corresponding to rDNA loci, identified in Navrátilová et al. 2022. (https://onlinelibrary.wiley.com/doi/full/10.1111/pbi.13816)
tc_df <- subset(tc_df, seqnames != "chr5H" & seqnames != "chr6H" | seqnames == "chr5H" & dominant_ctss<52608306 | seqnames == "chr5H" & dominant_ctss>53499223 | seqnames == "chr6H" & dominant_ctss<81918150 | seqnames == "chr6H" & dominant_ctss>82454047)

## ---- Differential expression analysis ---------------------------------------------------

ce$group <- factor(c(samples$group))
dds <- consensusClustersDESeq2(ce, ~group)
dds <- DESeq(dds)


