# 6.3 pre-exclusion DEA 6 times -------------------------------------------

# a version of 6.3 where I tried to run the code in parallel - seems slower if anything
# possibly because this computer doesn't really have much parallelt ability? 
# but this code could be useful if I want to run on Viking or something in future

# create matrix to loop through
# each row contains 24 random TCGA patient IDs from kansl1_wt
NTIMES_TO_RUN <- 6
NSAMPLES_TO_RUN <- 24
wt_mat <- matrix(sample(kansl1_wt, NSAMPLES_TO_RUN*NTIMES_TO_RUN), nrow = NTIMES_TO_RUN, ncol = NSAMPLES_TO_RUN)

# edits to run parallelised 
# convert matrix to list
wt_list <- split(wt_mat, row(wt_mat))
# register parallel background
register(SnowParam(workers = NTIMES_TO_RUN))

DESeq_ntimes <- function(wt, mutants, counts) {
  
  sample_ids <- c(unique(mutants), unique(wt))
  genotype <- c(rep("MUT", length(mutants)), 
                rep("WT", length(wt)))
  sample_info <- data.frame(row.names = sample_ids, genotype = genotype)
  
  # subset counts to only include the samples in sample info
  counts_subset <- counts[, sample_ids]
  counts_subset <- round(counts_subset)
  
  # create DESeq2 object
  dds <- DESeqDataSetFromMatrix(countData = counts_subset,
                                colData = sample_info,
                                design = ~ genotype)
  dds$genotype <- relevel(dds$genotype, ref = "WT")
  
  # run DESeq
  dds <- DESeq(dds)
  # store results
  dds_results <- results(dds)
  
  # change NA values to 0 and add a max adjusted p value
  dds_results$log2FoldChange[is.na(dds_results$log2FoldChange)] <- 0
  dds_results$padj[is.na(dds_results$padj) | dds_results$padj > 0.99] <- 0.99
  
  # add labels for colouring a volcano plot
  dds_results$DEA <- "NO" 
  dds_results$DEA[dds_results$log2FoldChange > 1 & dds_results$padj < 0.05] <- "UP"
  dds_results$DEA[dds_results$log2FoldChange < -1 & dds_results$padj < 0.05] <- "DOWN"
  
  # how many sig genes are there?
  sig_genes <- dds_results[which(dds_results$DEA %in% c("UP", "DOWN")), ]
  sig_genes <- rownames(sig_genes)
  return(length(sig_genes))
  
  
}

nsig24 <- 
  bplapply(X = wt_mat, FUN = DESeq_ntimes, 
           mutants = kansl1_muts, counts = counts)