# Capstone project analysis

# 1 Set-up ------------------------------------------------------------------

# load libraries
library(tidyverse)
library(dplyr)
library(maftools)
library(data.table)
library(googledrive)
library(ggplot2)
library(ggrepel)
library(DESeq2)
# library(fgsea) # couldn't be installed but ignoring this for now
library(BiocParallel)

# enable parallel computing
register(SerialParam())
# note: serialparam does not run parallelised cod
# it allows parallel code to be run unparallelised
# https://www.bioconductor.org/packages//release/bioc/vignettes/BiocParallel/inst/doc/Introduction_To_BiocParallel.html#quick-start

# log into google drive
googledrive::drive_auth()
# click 1

# 2 Identify tumours with KANSL1 mutation -----------------------------------

# read Whole Exome Sequencing maf
wxs_maf <- read.maf('data/WXS.maf')

# get summaries of all mutated genes
gene_tots <- getGeneSummary(wxs_maf)
gene_tots[grep('KANSL1', gene_tots$Hugo_Symbol),]$MutatedSamples
gene_tots[grep('KANSL1', gene_tots$Hugo_Symbol),]$AlteredSamples
# not sure what the difference between the mutated and altered samples columns are
# 24 patients had mutations in KANSL1

# lollipop plot could not be produced, cBioPortal will be used instead

# extract KANSL1 mutant patient IDs
kansl1_muts <- unique(wxs_maf@data[Hugo_Symbol == "KANSL1", Patient_Id])
kansl1_muts

# 3 Investigating consensus classifiers of KANSL1 mutants -------------------

# are the KANSL1 mutants more frequently in any one of the consensus classifiers?
con_class <- read.table("data/mRNA_gc47-TPMs_ConsensusClassifier.tsv", 
                        header = TRUE, sep = "\t")
summary(con_class)

# add Patient ID column with 01A removed from each of the IDs
con_class$Patient_ID <- substr(con_class$ID, 1, 12)

# what are the classifiers for the KANSL1 mutants we extracted?
con_class_kansl1 <- con_class[con_class$Patient_ID %in% kansl1_muts, ] |> 
  select('Patient_ID', 'consensusClass')

# 3.1 Create pie chart --------------------------------------------------------

# create two pie charts for visual comparison of whether or not there is a difference
# the proportion of classifiers in KANSL1 mutants and the entire cohort

png(file = 'plots/pies.png',
    width = 10, height = 6, units = 'in', res = 1000)
# set up plot matrix
par(mfrow = c(1,2))

# pie chart of entire cohort
class_counts_all <- table(con_class$consensusClass)
pie(class_counts_all,
    labels = names(class_counts_all),
    main = 'Consensus Classifiers of TCGA cohort')

# pie chart of KANSL1 mutants
class_counts_kansl1 <- table(con_class_kansl1$consensusClass)
pie(class_counts_kansl1,
    labels = names(class_counts_kansl1),
    main = 'Consensus Classifiers of KANSL1 mutants')

dev.off()

# most wedges are unchanged in size
# LumP wedge is a little larger in KANSL1 mutants, and Ba/Sq wedge a little smaller

# chi-squared to test if KANSL1 mutants are more likely to be LumP than expected

# add KANSL1 status to con_class
con_class$KANSL1 <- ifelse(con_class$Patient_ID %in% kansl1_muts, 'mutant', 'WT')
# add LumP status to con_class
con_class$LumP_status <- ifelse(con_class$consensusClass == 'LumP', 'TRUE', 'FALSE')

# create contingency table
table <- table(con_class$LumP_status, con_class$KANSL1)
table

# run chi-squared test
chisq.test(table)
# p-value 0.156
# no significant change 

# code could be adapted to investigate all 6
for (subtype in c('LumP', 'LumNS', 'Ba/Sq', 'Stroma-rich', 'NE-like', 'LumU')) { 
con_class$subtype_status <- ifelse(con_class$consensusClass == subtype, 'TRUE', 'FALSE')

# create contingency table
table <- table(con_class$subtype_status, con_class$KANSL1)

# run chi-squared test
print(chisq.test(table))
}
# none of them are significant
# also most of them are really small sample sizes so it wouldn't have mattered if they were

# tidy up con_class
con_class$LumP_status <- NULL
con_class$subtype_status <- NULL

# 4 DE analysis - KANSL1 mutants vs ALL non-mutants -----------------------------

# we need to exclude tumours who have mutations that might do the same thing as KANSL1 mutants
# these includes KANSL2 and WRD5 mutants (other subcomplexes of the NSL complex)
# as well as HAT1 and KDM1A

# extract patient IDs for all these genes
genes_to_extract <- c('KANSL2', 'WRD5', 'HAT1', 'KDM1A')

IDs_to_remove <- unique(wxs_maf@data[Hugo_Symbol %in% genes_to_extract, Patient_Id])
IDs_to_remove

# remove patient IDs who have mutations in these genes
kansl1_wt <- setdiff(kansl1_wt, IDs_to_remove)
length(kansl1_wt) # fewer than earlier - a success

# kansl1_muts is our object with KANSL1 mutant patient IDs
# create object of KANSL1 wildtype patient IDs
kansl1_wt <- setdiff(con_class$Patient_ID, kansl1_muts)
kansl1_wt

# create a sample info dataframe for differential expression
sample_ids <- c(kansl1_muts, kansl1_wt)
sample_ids

genotype <- c(rep("MUT", length(kansl1_muts)), 
              rep("WT", length(kansl1_wt)))
genotype

sample_info <- data.frame(row.names = sample_ids, genotype = genotype)

# read in mRNA data from online source (AM's google drive)
file_id <- '1djprP7DAcEwdUqugIOyQgM_JaYFmtcyl'
mrna_temp <- tempfile(fileext = ".tsv")
drive_download(as_id(file_id), path = mrna_temp, overwrite = TRUE)

# load RNAseq count data
counts <- read.table(mrna_temp, check.names = FALSE,
                     header = TRUE, row.names = 1, sep = '\t')

# remove temporary file
unlink(mrna_temp)
rm(mrna_temp)

# remove last 4 characters
colnames(counts) <- substr(colnames(counts), 1, nchar(colnames(counts)) - 4)

# sanity check - make sure samples in counts are the same as the samples taken from mutation data
counts <- counts[ , sample_ids]

# round the counts and only look at complete cases (i.e. all genes, I think?)
counts <- round(counts[complete.cases(counts), ])

# create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = sample_info,
                              design = ~ genotype)

dds$genotype <- relevel(dds$genotype, ref = "WT")
dds <- DESeq(dds)

# store results
dds_results <- results(dds)
dds_results

# change NA values to 0 and add a max adjusted p value
dds_results$log2FoldChange[is.na(dds_results$log2FoldChange)] <- 0
dds_results$padj[is.na(dds_results$padj) | dds_results$padj > 0.99] <- 0

# add labels for colouring a volcano plot
dds_results$DEA <- "NO" 
dds_results$DEA[dds_results$log2FoldChange > 1 & dds_results$padj < 0.05] <- "UP"
dds_results$DEA[dds_results$log2FoldChange < -1 & dds_results$padj < 0.05] <- "DOWN"

# add gene symbols as column for easy plot labelling
dds_results$symbol <- rownames(dds_results)

# create pi values (fold change multiplied by stat significance) for later gene ranking
dds_results$pi <- dds_results$log2FoldChange * -log10(dds_results$padj)
dds_results_pi_sorted <- dds_results[order(dds_results$pi),]

# create reduced dataframe of most significantly different genes, for labelling
# using pi values means you prioritise most biologically significant
top_genes <- c(head(dds_results_pi_sorted$symbol, 20), tail(dds_results_pi_sorted$symbol, 20))
genes_to_label <- dds_results[dds_results$symbol %in% top_genes, ]

# create a labelled, coloured and annotated volcano plot
png(file = 'plots/DEA_KANSL1_wt.png',
    width = 10, height = 6, units = 'in', res = 1000)

ggplot(dds_results, aes(x=log2FoldChange, y=-log10(padj))) + 
  geom_point(aes(colour = DEA), show.legend = FALSE) + 
  scale_colour_manual(values = c("blue", "gray", "red")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted") +
  geom_vline(xintercept = c(-1,1), linetype = "dotted") + 
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  geom_text_repel(size=2, data=genes_to_label, aes(x=log2FoldChange, y=-log10(padj), label=symbol), max.overlaps = Inf)

dev.off()

# 5 DE analysis - 24vs24 ---------------------------------------------------

# kansl1_muts is our object with KANSL mutant patient IDs
# create object of KANSL1 wildtype patient IDs
kansl1_wt <- setdiff(con_class$Patient_ID, kansl1_muts)
kansl1_wt

# we need to exclude tumours who have mutations that might do the same thing as KANSL1 mutants
# these includes KANSL2 and WRD5 mutants (other subcomplexes of the NSL complex)
# as well as HAT1 and KDM1A

# extract patient IDs for all these genes
genes_to_extract <- c('KANSL2', 'WRD5', 'HAT1', 'KDM1A')

IDs_to_remove <- unique(wxs_maf@data[Hugo_Symbol %in% genes_to_extract, Patient_Id])
IDs_to_remove

# remove patient IDs who have mutations in these genes
kansl1_wt <- setdiff(kansl1_wt, IDs_to_remove)
length(kansl1_wt) # fewer than earlier - a success

# create object of 24 randomly selected KANSL1 wildtype patient IDs
kansl1_wt <- sample(kansl1_wt, 24)

# create a sample info dataframe for differential expression
sample_ids <- c(kansl1_muts, kansl1_wt)
sample_ids

genotype <- c(rep("MUT", length(kansl1_muts)), 
              rep("WT", length(kansl1_wt)))
genotype

sample_info <- data.frame(row.names = sample_ids, genotype = genotype)

# read in mRNA data from online source (AM's google drive)
file_id <- '1djprP7DAcEwdUqugIOyQgM_JaYFmtcyl'
mrna_temp <- tempfile(fileext = ".tsv")
drive_download(as_id(file_id), path = mrna_temp, overwrite = TRUE)

# load RNAseq count data
counts <- read.table(mrna_temp, check.names = FALSE,
                     header = TRUE, row.names = 1, sep = '\t')

# remove temporary file
unlink(mrna_temp)
rm(mrna_temp)

# adapt column names (remove last four characters)
colnames(counts) <- substr(colnames(counts), 1, nchar(colnames(counts)) - 4)

# sanity check - make sure samples in counts are the same as the samples taken from mutation data
counts <- counts[ , sample_ids]

# round the counts and only look at complete cases (i.e. all genes, I think?)
counts <- round(counts[complete.cases(counts), ])

# create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = sample_info,
                              design = ~ genotype)

dds$genotype <- relevel(dds$genotype, ref = "WT")
dds <- DESeq(dds)

# store results
dds_results <- results(dds)
dds_results

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
length(sig_genes)
# 

# add gene symbols as column for easy plot labelling
dds_results$symbol <- rownames(dds_results)

# create pi values (fold change multiplied by stat significance) for later gene ranking
dds_results$pi <- dds_results$log2FoldChange * -log10(dds_results$padj)
dds_results_pi_sorted <- dds_results[order(dds_results$pi),]

# create reduced dataframe of most significantly different genes, for labelling
# using pi values means you prioritise most biologically significant
top_genes <- c(head(dds_results_pi_sorted$symbol, 20), tail(dds_results_pi_sorted$symbol, 20))
genes_to_label <- dds_results[dds_results$symbol %in% top_genes, ]

# create volcano plot
png(file = 'plots/DEA_KANSL1_wt24.png',
    width = 8, height = 5, units = 'in', res = 1000)

# create a labelled, coloured and annotated volcano plot
ggplot(dds_results, aes(x=log2FoldChange, y=-log10(padj))) + 
  geom_point(aes(colour = DEA),
             show.legend = FALSE) + 
  scale_colour_manual(values = c("blue", "gray", "red")) +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dotted") +
  geom_vline(xintercept = c(-1,1),
             linetype = "dotted") + 
  theme_classic() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")) +
  geom_text_repel(size=2,
                  data=genes_to_label,
                  aes(x=log2FoldChange, y=-log10(padj),label=symbol),
                  max.overlaps = Inf)

dev.off()

# 5.1 Savepoint A -------------------------------------------------------------

# create save point so you can come back to the same random 24
# save.image("savepointA_DEA1-complete.RData")

# load("savepointA_DEA1-complete.RData")

# 6 DEA excluding some mutants ------------------------------------------------

# create data frame which PolyPhen and SIFT scores can be added to (taken from GenomeNexus)
mut_info <- unique(wxs_maf@data[Hugo_Symbol == "KANSL1", Patient_Id, Protein_Change])
mut_info <- mut_info[order(mut_info$Patient_Id), ] # order Patient IDs alphabetically
mut_info$PolyPhen <- c('NA', 'NA', 'NA', 'NA', 0.03, 'NA','NA', 'NA', 0.99, 0.99, 'NA', 'NA', 0.24, 'NA', 'NA', 0.73, 'NA', 0.99, 0.91, 'NA', 0.66, 1.00, 0.12, 0.04, 'NA')
mut_info$SIFT <- c('NA', 'NA', 'NA', 'NA', 0.18, 'NA', 'NA', 'NA', 0.00, 0.00, 'NA', 'NA', 0.16, 'NA', 'NA', 0.01, 'NA', 0.14, 0.01, 'NA', 0.01, 0.00, 0.10, 0.26, 'NA')

# filter KANSL1 mutations using thresholds consistent with those used by Ensembl
mut_info <- mut_info |> filter(PolyPhen == 'NA' | PolyPhen > 0.446 | SIFT == 'NA' | SIFT < 0.05)
nrow(mut_info)
# 21 mutants now remain

# 6.1 Does excluding genes make a difference? ---------------------------------

# excluding mutants which don't do anything should reduce increase the number of
# significantly differentially expressed genes

# it seems that the number of significantly differentially expressed genes can vary
# quite a lot, so I will run differential expression before and after gene exclusion 5 times

# 6.2 Set up to loop DEA --------------------------------------------------

# this is much of the same set up as earlier

# get KANSL1 mutants again
kansl1_muts <- unique(wxs_maf@data[Hugo_Symbol == "KANSL1", Patient_Id])
kansl1_muts

# kansl1_muts is our object with KANSL mutant patient IDs
# create object of KANSL1 wildtype patient IDs
kansl1_wt <- setdiff(con_class$Patient_ID, kansl1_muts)
kansl1_wt

# we need to exclude tumours who have mutations that might do the same thing as KANSL1 mutants
# these includes KANSL2 and WRD5 mutants (other subcomplexes of the NSL complex)
# as well as HAT1 and KDM1A

# extract patient IDs for all these genes
genes_to_extract <- c('KANSL2', 'WRD5', 'HAT1', 'KDM1A')

IDs_to_remove <- unique(wxs_maf@data[Hugo_Symbol %in% genes_to_extract, Patient_Id])
IDs_to_remove

# remove patient IDs who have mutations in these genes
kansl1_wt <- setdiff(kansl1_wt, IDs_to_remove)
length(kansl1_wt) # fewer than earlier - a success

# read in mRNA data from online source (AM's google drive)
file_id <- '1djprP7DAcEwdUqugIOyQgM_JaYFmtcyl'
mrna_temp <- tempfile(fileext = ".tsv")
drive_download(as_id(file_id), path = mrna_temp, overwrite = TRUE)

# load RNAseq count data
counts <- read.table(mrna_temp, check.names = FALSE,
                     header = TRUE, row.names = 1, sep = '\t')

# remove temporary file
unlink(mrna_temp)
rm(mrna_temp)

# adapt column names (remove last four characters)
colnames(counts) <- substr(colnames(counts), 1, nchar(colnames(counts)) - 4)

# round the counts and only look at complete cases (i.e. all genes, I think?)
counts <- counts[complete.cases(counts), ]


# 6.3 pre-exclusion DEA 6 times -------------------------------------------

# create matrix to loop through
# each row contains 24 random TCGA patient IDs from kansl1_wt
NTIMES_TO_RUN <- 12
NSAMPLES_TO_RUN <- 24
wt_mat <- matrix(sample(kansl1_wt, NSAMPLES_TO_RUN*NTIMES_TO_RUN),
                 nrow = NTIMES_TO_RUN, ncol = NSAMPLES_TO_RUN)

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

nsig24 <- apply(X = wt_mat, MARGIN = 1, FUN = DESeq_ntimes, 
        mutants = kansl1_muts, counts = counts)

saveRDS(nsig24, file = 'outputs/nsig24.RDS')

# 6.4 DEA with some mutants excluded, 6 times -----------------------------------------

# resetting KANSL1 mutant object - won't include mutants which don't break KANSL1
kansl1_muts <- unique(mut_info$Patient_Id)
kansl1_muts

# create wildtype KANSL1 object
kansl1_wt <- setdiff(con_class$Patient_ID, kansl1_muts)
length(kansl1_wt)

# remove some wt mutants, as done in 6.2
kansl1_wt <- setdiff(kansl1_wt, IDs_to_remove)
length(kansl1_wt)

# create matrix of KANSL1 wildtypes to loop through
NTIMES_TO_RUN <- 12
NSAMPLES_TO_RUN <- 20 # kansl1_muts has gone down to 20 so I have set this to do the same? I guess because I added a unique somewhere?
wt_mat <- matrix(sample(kansl1_wt, NSAMPLES_TO_RUN*NTIMES_TO_RUN),
                 nrow = NTIMES_TO_RUN, ncol = NSAMPLES_TO_RUN)

# run DESeq n times and save
nsig20 <- apply(X = wt_mat, MARGIN = 1, FUN = DESeq_ntimes, 
                mutants = kansl1_muts, counts = counts)

saveRDS(nsig20, file = 'outputs/nsig20.RDS')

# 6.5 Looking at the difference in nsig genes -----------------------------------


# create tidy data frame for plotting a box plot
nsig_all <- data.frame(nsig24)
nsig_all$KANSL1_mutants <- 'All'
nsig_all$nsig_genes <- nsig_all$nsig24
nsig_all$nsig24 <- NULL
nsig_all

nsig_excl <- data.frame(nsig20)
nsig_excl$KANSL1_mutants <- 'Excluded'
nsig_excl$nsig_genes <- nsig_excl$nsig20
nsig_excl$nsig20 <- NULL
nsig_excl

nsigs <- rbind(nsig_all, nsig_excl)

saveRDS(nsigs, file = 'outputs/nsigs.RDS')

# some helpful summary stastics
nsigs_summary <- nsigs |> 
  group_by(KANSL1_mutants) |> 
  summarise(mean = mean(nsig_genes),
            median = median(nsig_genes),
            sd = sd(nsig_genes),
            n = length(nsig_genes),
            se = sd/sqrt(n))

# create plot
png(file = 'plots/excluding_mutants.png',
    width = 5, height = 7, units = 'in', res = 1000)

ggplot() +
  geom_point(data = nsigs, aes(x = KANSL1_mutants, y = nsig_genes),
             position = position_jitter(width = 0.1, height = 0)) + 
  geom_errorbar(data = nsigs_summary,
                aes(x = KANSL1_mutants, ymin = mean - se, ymax = mean + se),
                width = 0.5) +
  geom_errorbar(data = nsigs_summary,
              aes(x = KANSL1_mutants, ymin = mean, ymax = mean),
              width = 0.3) +
  annotate("segment", x = 1, xend = 2, y = 1300, yend = 1300,
           colour = "black") +
  annotate('text', x = 1.5, y = 1350,
           label = 'n.s.') +
  theme_classic()

dev.off()

# looks like the reverse of what was expected

# run t test
t.test(formula = nsig_genes ~ KANSL1_mutants, data = nsigs)

# 6.6 Savepoint B ---------------------------------------------------------

# save.image("savepointB_DEA2-complete.RData")


# 7 New semester: DEA KANSL1 mutants vs ALL KANSL1 WT cancers ---------------

# BUT this time performing some sensible exclusions first
# Read Lab Book to understand the reasoning behind the new strategy
# new year new DEA

# 7.1 exclusions ----------------------------------------------------------

# exclude tumours with mutations that might do the same thing as KANSL1 mutants
# KANSL2, WDR5, HAT1 and KDM1A

# extract patient IDs for all these genes
genes_to_extract <- c('KANSL2', 'WRD5', 'HAT1', 'KDM1A')

IDs_to_remove <- unique(wxs_maf@data[Hugo_Symbol %in% genes_to_extract, Patient_Id])
IDs_to_remove

# kansl1_muts is our object with KANSL1 mutant patient IDs
# create object of KANSL1 wildtype patient IDs
kansl1_wt <- setdiff(con_class$Patient_ID, kansl1_muts)
length(kansl1_wt) # 384

# remove patient IDs who have mutations in these genes
kansl1_wt <- setdiff(kansl1_wt, IDs_to_remove)
length(kansl1_wt) # 374 - fewer than earlier - a success

# exclude KANSL1 mutants where the mutation might not actually do anything
# keep high PolyPhen values and and low SIFT values

# create data frame which PolyPhen and SIFT scores can be added to (taken from GenomeNexus)
mut_info <- unique(wxs_maf@data[Hugo_Symbol == "KANSL1", Patient_Id, Protein_Change])
mut_info <- mut_info[order(mut_info$Patient_Id), ] # order Patient IDs alphabetically
mut_info$PolyPhen <- c('NA', 'NA', 'NA', 'NA', 0.03, 'NA','NA', 'NA', 0.99, 0.99, 'NA', 'NA', 0.24, 'NA', 'NA', 0.73, 'NA', 0.99, 0.91, 'NA', 0.66, 1.00, 0.12, 0.04, 'NA')
mut_info$SIFT <- c('NA', 'NA', 'NA', 'NA', 0.18, 'NA', 'NA', 'NA', 0.00, 0.00, 'NA', 'NA', 0.16, 'NA', 'NA', 0.01, 'NA', 0.14, 0.01, 'NA', 0.01, 0.00, 0.10, 0.26, 'NA')

# filter KANSL1 mutations using thresholds consistent with those used by Ensembl
mut_info <- mut_info |> filter(PolyPhen == 'NA' | PolyPhen > 0.446 | SIFT == 'NA' | SIFT < 0.05)
nrow(mut_info)
# there are now 21 KANSL1 mutants remaining

# setting KANSL1 mutant object
kansl1_muts <- unique(mut_info$Patient_Id)
kansl1_muts
length(kansl1_muts)

# clear up
rm(genes_to_extract, mut_info, IDs_to_remove)

# 7.2 DEA -----------------------------------------------------------------

# create a sample info dataframe for differential expression
sample_ids <- c(kansl1_muts, kansl1_wt)
sample_ids

genotype <- c(rep("MUT", length(kansl1_muts)), 
              rep("WT", length(kansl1_wt)))
genotype

sample_info <- data.frame(row.names = sample_ids, genotype = genotype)

# clear up
rm(genotype)

# read in mRNA data from online source (AM's google drive)
file_id <- '1djprP7DAcEwdUqugIOyQgM_JaYFmtcyl'
mrna_temp <- tempfile(fileext = ".tsv")
drive_download(as_id(file_id), path = mrna_temp, overwrite = TRUE)

# load RNAseq count data
counts <- read.table(mrna_temp, check.names = FALSE,
                     header = TRUE, row.names = 1, sep = '\t')

# remove temporary file
unlink(mrna_temp)
rm(mrna_temp, file_id)

# remove last 4 characters
colnames(counts) <- substr(colnames(counts), 1, nchar(colnames(counts)) - 4)

# sanity check - make sure samples in counts are the same as the samples taken from mutation data
counts <- counts[ , sample_ids]

# round the counts and only look at complete cases (i.e. all genes, I think?)
counts <- round(counts[complete.cases(counts), ])

# create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = sample_info,
                              design = ~ genotype)

dds$genotype <- relevel(dds$genotype, ref = "WT")
dds <- DESeq(dds)

# store results
dds_results <- results(dds)
dds_results

# tidy up
rm(dds)

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
length(sig_genes)
# 2328 

# add gene symbols as column for easy plot labelling
dds_results$symbol <- rownames(dds_results)

# create pi values (fold change multiplied by stat significance) for later gene ranking
dds_results$pi <- dds_results$log2FoldChange * -log10(dds_results$padj)
dds_results_pi_sorted <- dds_results[order(dds_results$pi),]

# create reduced dataframe of most significantly different genes, for labelling
# using pi values means you prioritise most biologically significant
top_genes <- c(head(dds_results_pi_sorted$symbol, 20), tail(dds_results_pi_sorted$symbol, 20))
genes_to_label <- dds_results[dds_results$symbol %in% top_genes, ]

# create volcano plot
png(file = 'plots/DEA_kansl1.png',
    width = 8, height = 5, units = 'in', res = 1000)

# create a labelled, coloured and annotated volcano plot
ggplot(dds_results, aes(x=log2FoldChange, y=-log10(padj))) + 
  geom_point(aes(colour = DEA),
             show.legend = FALSE) + 
  scale_colour_manual(values = c("blue", "gray", "red")) +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dotted") +
  geom_vline(xintercept = c(-1,1),
             linetype = "dotted") + 
  theme_classic() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")) +
  geom_text_repel(size=2,
                  data=genes_to_label,
                  aes(x=log2FoldChange, y=-log10(padj),label=symbol),
                  max.overlaps = Inf)

dev.off()
