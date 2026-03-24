# Capstone project analysis

# 1 Set-up ------------------------------------------------------------------

# load packages
library(tidyverse)
library(dplyr)
library(maftools)
library(data.table)
library(googledrive)
library(ggplot2)
library(ggrepel)
library(DESeq2)
library(fgsea)
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
# run up to here for set up of later analysis

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
# add some more info for a summary table too
for (subtype in c('LumP', 'LumNS', 'Ba/Sq', 'Stroma-rich', 'NE-like', 'LumU')) { 
con_class$subtype_status <- ifelse(con_class$consensusClass == subtype, 'TRUE', 'FALSE')

# create contingency table
table <- table(con_class$subtype_status, con_class$KANSL1)

# run chi-squared test
print(subtype)
print(chisq.test(table))

# get some more statistics for a summary table
# the number of each subtype in the TCGA cohort
# the number of each subtype in KANSL1 mutants
# the proportion of each subtype in the TCGA cohort
# the proportion of each subtype in KANSL1 mutants
no_TCGA <- sum(con_class$subtype_status == 'TRUE')
no_KANSL1muts <- sum(con_class$subtype_status == 'TRUE' & con_class$KANSL1 == 'mutant')
prop_TCGA <- no_TCGA/nrow(con_class)
prop_KANSL1muts <- no_KANSL1muts/sum(con_class$KANSL1 == 'mutant')
print(no_TCGA)
print(no_KANSL1muts)
print(prop_TCGA)
print(prop_KANSL1muts)
  
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
# 21 notable mutations remaining

# setting KANSL1 mutant object
kansl1_muts <- unique(mut_info$Patient_Id)
kansl1_muts
length(kansl1_muts)
# there's 20 KANSL1 mutants left - two of the mutations identified above are from the same patient

# clear up
rm(genes_to_extract, mut_info, IDs_to_remove)

# 7.2 DEA set up -----------------------------------------------------------------

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

# tidy up
rm(mrna_temp, file_id)

# remove last 4 characters
colnames(counts) <- substr(colnames(counts), 1, nchar(colnames(counts)) - 4)

# sanity check - make sure samples in counts are the same as the samples taken from mutation data
counts <- counts[ , sample_ids]

# round the counts and only look at complete cases (i.e. all genes, I think?)
counts <- round(counts[complete.cases(counts), ])

# 7.3 run DEA -----------------------------------------------------------------

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
# 459 

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
# png(file = 'plots/DEA_kansl1.png',
#    width = 8, height = 5, units = 'in', res = 1000)
# hidden to avoid accidental edits

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

# dev.off()

# save data
dds_df <- as.data.frame(dds_results)
write.csv(dds_df, "results/KANSL1_mutvsWT_DEA_results.csv", row.names = FALSE, col.names = TRUE)


# 7.3.1 Is KANSL1 differentially expressed--------------------------------------------------

# this feels like it might be a prudent thing to check - that KANSL1 is actually 
# down in KANSL1 mutants. also helps us to orientate ourselves

# get KANSL1 DEA info
#dds_results[dds_results$symbol == "KANSL1", ]
# hm barely no change
# should this be a concern?

# volcano plot with KANSL1 labelled
#genes_to_label <- dds_results[dds_results$symbol %in% c(top_genes, 'KANSL1'), ]

#png(file = 'plots/DEA_kansl1-labelled.png',
 #   width = 8, height = 5, units = 'in', res = 1000)

#ggplot(dds_results, aes(x=log2FoldChange, y=-log10(padj))) + 
 # geom_point(aes(colour = DEA),
  #           show.legend = FALSE) + 
  #scale_colour_manual(values = c("blue", "gray", "red")) +
  #geom_hline(yintercept = -log10(0.05),
  #           linetype = "dotted") +
  #geom_vline(xintercept = c(-1,1),
  #           linetype = "dotted") + 
  #theme_classic() +
  #theme(panel.border = element_blank(),
  #      panel.grid.major = element_blank(),
  #      panel.grid.minor = element_blank(),
  #      axis.line = element_line(colour = "black")) +
  #geom_text_repel(size=2,
  #                data=genes_to_label,
  #                aes(x=log2FoldChange, y=-log10(padj),label=symbol),
  #                max.overlaps = Inf)

# dev.off()

# good news - this isn't important to check
# hiding this code with hashes


# 7.3.2 Checking which way round the volcano plot is -------------------------------------------------------------------

# this is alos a prudent thing to check
# pick a gene from the volcano plot, look at the TPMs in counts
counts['MOG', ]
# mostly 0, 1 and 2 TPMs
# some much higher - TCGA-XF-AAMG has 95
# is this a KANSL1 mutant?
sample_info['TCGA-XF-AAMG', ]
# YES
# so feeling pretty good that the red side is genes up in KANSL1 mutants
# can look at all KANSL1 mutants pretty easily
counts['MOG' , kansl1_muts]
# definitely higher than background 0s and 1s

# can do the same with a blue gene
counts['TENM2', ]
counts['TENM2', kansl1_muts]
# MOST of these are just in double digits, where the others as a whole were in triple/quadruple

# SO I'm pretty happy - reds are up in KANSL1 mutants, blues are down in KANSL1 mutants
# and we definitely set ref = WT, which backs this up

# 7.4 GSEA on KANSL1 mutant vs WT DEA --------------------------------------------

# load MSigDB gene set
genesets = gmtPathways("data/c2.all.v2023.2.Hs.symbols.gmt")

# use pi values to rank genes from "most up" to "most down" in the comparison
prerank <- dds_results[c("symbol", "pi")]
prerank <- setNames(prerank$pi, prerank$symbol)
str(prerank)

# run fgsea
fgseaRes <- fgsea(pathways = genesets, stats = prerank, minSize=15, maxSize = 500)

# store top10 most enriched
# positive enrichment is (relatively) up in the KANSL1 mutant tumours

#~ negative enrichment is (relatively) up in the KANSL1 wildtype tumours
top10_fgseaRes <- head(fgseaRes[order(pval), ], 10)
top10_fgseaRes
# top 7 have padj <0.05

# these seem to keep changing each time the code is run
# look at more - 100 
top100_fgseaRes <- head(fgseaRes[order(pval), ], 100)
top100_fgseaRes

# there's a lot tied - let's just look at those with a padj <0.05 for now
sig_fgseaRes <- fgseaRes[padj < 0.05, ]
  
# create bar chart of normalised enrichment scores (NES) for top10 hits
# hidden to avoid overwriting
# png(file = 'plots/FGSEA_KANSL1_mut_vs_WT_padj0.05.png',
#    width = 12, height = 5, units = 'in', res = 1000)

ggplot(sig_fgseaRes, aes(x = NES, y=reorder(pathway, -pval), fill = factor(sign(NES)))) + 
  geom_bar(stat = "identity", width = 0.8) +
  labs(title = "GSEA", x = "Normalised Enrichment Score (NES)", y = "Pathway") +
  theme_minimal(base_size = 16) +
  scale_fill_manual(values = c("#0754A2", "#B10029"), guide = "none") +
  scale_y_discrete(labels = function(x) gsub("^HALLMARK_", "", x)) +
  theme(axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# dev.off()

# 7.4.1 Looking at Keratinisation -------------------------------------------------------------------

# what are the leading edge genes for the keratinisation gene set?
krt_lEdge <- unlist(sig_fgseaRes[pathway == 'REACTOME_KERATINIZATION', leadingEdge])

# replot volcano plot with these genes labelled to see where they are
genes_to_label <- dds_results[dds_results$symbol %in% krt_lEdge, ]

# png(file = 'plots/DEA_keratins-labelled.png',
#   width = 8, height = 5, units = 'in', res = 1000)

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

# dev.off()

# 8 Running 20vs20 DEA ----------------------------------------------------

# it would back up any interesting genes we've found to run the DEA 20vs20
# (20 KANSL1 mutants vs 20 WT). since the 20 WTs are random and the genes that come up
# are pretty different each time, I want to write some code which runs the DEA multiple times,
# stores the results from each run, and calculates a the proportion of times a gene was significant

# set up matrix of WT selection to run
NTIMES_TO_RUN <- 3
NSAMPLES_TO_RUN <- 20
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
  
  # add UP/DOWN regulated labels
  dds_results$DEA <- "NO" 
  dds_results$DEA[dds_results$log2FoldChange > 1 & dds_results$padj < 0.05] <- "UP"
  dds_results$DEA[dds_results$log2FoldChange < -1 & dds_results$padj < 0.05] <- "DOWN"
  
  # return results
  return(as.data.frame(dds_results))
  
}

# create emptor vector (type = list) to store results in 
DEA_results <- vector("list", NTIMES_TO_RUN)

for (i in seq_len(NTIMES_TO_RUN)) {
  
  # select WT samples for this run
  wt_samples <- wt_mat[i, ]
  
  # run DESeq for this run and store run results in res
  res <- DESeq_ntimes(wt = wt_samples, mutants = kansl1_muts, counts = counts)
  
  # add gene name in a column which isn't the row name
  # done here to avoid naming complications
  res$gene <- rownames(res)
  
  # add column saying which run this is
  res$run <- i
  
  # add run results to overall results vector
  DEA_results[[i]] <- res
  
}

# combine the DEA results by row (gene)
DEA_combined <- do.call(rbind, DEA_results)

# tidy up - DEA_results is huge
# rm(DEA_results)

# create summary table
DEA_summary <- DEA_combined |> 
  group_by(gene) |> 
  summarise(n_UP = sum(DEA == "UP"),
    n_DOWN = sum(DEA == "DOWN"),
    prop_UP = n_UP / NTIMES_TO_RUN,
    prop_DOWN = n_DOWN / NTIMES_TO_RUN,
    prop_DE = (n_UP + n_DOWN) / NTIMES_TO_RUN,
    mean_log2FC = mean(log2FoldChange),
    median_log2FC = median(log2FoldChange),
    sd_log2FC = sd(log2FoldChange))

# how does this look
head(DEA_summary, n = 10)

# extract reproducible genes which are significant > half the time
rep_genes <- DEA_summary[DEA_summary$prop_DE > 0.5, ]
dim(rep_genes)

# this code is now ready for a larger number of runs another time


# 9 Investigating keratinisation ----------------------------------

# we'd expect to see more keratinisation in Ba/Sq subclass
# might KANSL1 have any sort of role in regulating this?

# 9.1 DEA Ba/Sq mutants vs other mutants -------------------------------------------------------

# can run from here instead of 7

# exclude mutants that prob don't do anything
mut_info <- unique(wxs_maf@data[Hugo_Symbol == "KANSL1", Patient_Id, Protein_Change])
mut_info <- mut_info[order(mut_info$Patient_Id), ] # order Patient IDs alphabetically
mut_info$PolyPhen <- c('NA', 'NA', 'NA', 'NA', 0.03, 'NA','NA', 'NA', 0.99, 0.99, 'NA', 'NA', 0.24, 'NA', 'NA', 0.73, 'NA', 0.99, 0.91, 'NA', 0.66, 1.00, 0.12, 0.04, 'NA')
mut_info$SIFT <- c('NA', 'NA', 'NA', 'NA', 0.18, 'NA', 'NA', 'NA', 0.00, 0.00, 'NA', 'NA', 0.16, 'NA', 'NA', 0.01, 'NA', 0.14, 0.01, 'NA', 0.01, 0.00, 0.10, 0.26, 'NA')
mut_info <- mut_info |> filter(PolyPhen == 'NA' | PolyPhen > 0.446 | SIFT == 'NA' | SIFT < 0.05)
nrow(mut_info)
kansl1_muts <- unique(mut_info$Patient_Id)
kansl1_muts
rm(mut_info)

# create vector of Ba/Sq patient IDs with kansl1 mutations, and one with all other kansl1 mutants
basq_patients <- con_class$Patient_ID[con_class$consensusClass == 'Ba/Sq']
basq_kansl1 <- intersect(kansl1_muts, basq_patients)
other_kansl1 <- setdiff(kansl1_muts, basq_kansl1)

# run DEA, with other KANSL1 mutatns set as reference
# create a sample info dataframe for differential expression
sample_ids <- c(basq_kansl1, other_kansl1)
class <- c(rep("BaSq", length(basq_kansl1)), 
              rep("other", length(other_kansl1)))
sample_info <- data.frame(row.names = sample_ids, class = class)
rm(class)
sample_info

# read in mRNA data
file_id <- '1djprP7DAcEwdUqugIOyQgM_JaYFmtcyl'
mrna_temp <- tempfile(fileext = ".tsv")
drive_download(as_id(file_id), path = mrna_temp, overwrite = TRUE)
counts <- read.table(mrna_temp, check.names = FALSE,
                     header = TRUE, row.names = 1, sep = '\t')
unlink(mrna_temp)
rm(mrna_temp, file_id)

# sort counts object
colnames(counts) <- substr(colnames(counts), 1, nchar(colnames(counts)) - 4)
counts <- counts[ , sample_ids]
counts <- round(counts[complete.cases(counts), ])

# run DEseq
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = sample_info,
                              design = ~ class)

dds$class <- relevel(dds$class, ref = "other")
dds <- DESeq(dds)
dds_results <- results(dds)
dds_results
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
# 1401 - that's lots! a good thing? 

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
# png(file = 'plots/DEA_BaSq_KANSL1_muts_vs_others.png',
#    width = 8, height = 5, units = 'in', res = 1000)
# hidden to avoid accidental edits

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

# dev.off()

# lots of keratins are up in the Ba/Sq HOWEVER this doesn't yet show much
# we know that Ba/Sq have lots of keratins regardless so not necessarily a role for KANSL1 here

# 9.2 DEA no Ba/Sq mutants or WTs -------------------------------------------------------

# can run from here instead of 7

# 9.2.1 Exclusions --------------------------------------------------------

# this is to make sure the changes in keratinisation aren't just because of the fewer KANSL1 mutants 
# which are Ba/Sq, as observed in the pie charts

# create KANSL1 WT set of IDs - make same exclusions as earlier
genes_to_extract <- c('KANSL2', 'WRD5', 'HAT1', 'KDM1A')
IDs_to_remove <- unique(wxs_maf@data[Hugo_Symbol %in% genes_to_extract, Patient_Id])
kansl1_wt <- setdiff(con_class$Patient_ID, kansl1_muts)
length(kansl1_wt) # 384
# remove patient IDs who have mutations in these genes
kansl1_wt <- setdiff(kansl1_wt, IDs_to_remove)
length(kansl1_wt) # 374
# now remove patient IDs for WT Ba/Sq tumours
basq_patients <- con_class$Patient_ID[con_class$consensusClass == 'Ba/Sq']
basq_kansl1_wt <- intersect(kansl1_wt, basq_patients)
kansl1_wt <- setdiff(kansl1_wt, basq_kansl1_wt)
length(kansl1_wt)
rm(IDs_to_remove, genes_to_extract)

# create KANSL1 mutant set of IDs - make same exclusions as earlier
mut_info <- unique(wxs_maf@data[Hugo_Symbol == "KANSL1", Patient_Id, Protein_Change])
mut_info <- mut_info[order(mut_info$Patient_Id), ] # order Patient IDs alphabetically
mut_info$PolyPhen <- c('NA', 'NA', 'NA', 'NA', 0.03, 'NA','NA', 'NA', 0.99, 0.99, 'NA', 'NA', 0.24, 'NA', 'NA', 0.73, 'NA', 0.99, 0.91, 'NA', 0.66, 1.00, 0.12, 0.04, 'NA')
mut_info$SIFT <- c('NA', 'NA', 'NA', 'NA', 0.18, 'NA', 'NA', 'NA', 0.00, 0.00, 'NA', 'NA', 0.16, 'NA', 'NA', 0.01, 'NA', 0.14, 0.01, 'NA', 0.01, 0.00, 0.10, 0.26, 'NA')
mut_info <- mut_info |> filter(PolyPhen == 'NA' | PolyPhen > 0.446 | SIFT == 'NA' | SIFT < 0.05)
nrow(mut_info)
kansl1_muts <- unique(mut_info$Patient_Id)
kansl1_muts
rm(mut_info)

# now exclude KANSL1 mutants which are Ba/Sq
basq_kansl1_muts <- intersect(kansl1_muts, basq_patients)
length(basq_kansl1_muts)
kansl1_muts <- setdiff(kansl1_muts, basq_kansl1_muts)
length(kansl1_muts)
# tidy up
rm(basq_kansl1_muts, basq_kansl1_wt, basq_patients)

# 9.2.2 Set Up DEA --------------------------------------------------------

# create a sample info dataframe for differential expression
sample_ids <- c(kansl1_muts, kansl1_wt)
genotype <- c(rep("MUT", length(kansl1_muts)), 
              rep("WT", length(kansl1_wt)))
sample_info <- data.frame(row.names = sample_ids, genotype = genotype)
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
# tidy up
rm(mrna_temp, file_id)

# transform counts - remove last 4 characters
colnames(counts) <- substr(colnames(counts), 1, nchar(colnames(counts)) - 4)
# sanity check - make sure samples in counts are the same as the samples taken from mutation data
counts <- counts[ , sample_ids]
# round the counts and only look at complete cases (i.e. all genes, I think?)
counts <- round(counts[complete.cases(counts), ])

# 9.2.3 Run DEA -----------------------------------------------------------

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
length(sig_genes)# 459 

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
# png(file = 'plots/DEA_kansl1_no_BaSq.png', width = 8, height = 5, units = 'in', res = 1000)
# hidden to avoid accidental edits

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

# dev.off()

# save data
dds_df <- as.data.frame(dds_results)
# write.csv(dds_df, "results/KANSL1_mutvsWT_DEA_results_excl_BaSq.csv", row.names = FALSE, col.names = TRUE)

# 9.2.4 GSEA --------------------------------------------------------------

# load MSigDB gene set
genesets = gmtPathways("data/c2.all.v2023.2.Hs.symbols.gmt")

# use pi values to rank genes from "most up" to "most down" in the comparison
prerank <- dds_results[c("symbol", "pi")]
prerank <- setNames(prerank$pi, prerank$symbol)
str(prerank)

# run fgsea
fgseaRes <- fgsea(pathways = genesets, stats = prerank, minSize=15, maxSize = 500)

# store top10 most enriched
# positive enrichment is (relatively) up in the KANSL1 mutant tumours
# negative enrichment is (relatively) up in the KANSL1 wildtype tumours
top10_fgseaRes <- head(fgseaRes[order(pval), ], 10)
top10_fgseaRes

# isolate those with a padj <0.05
sig_fgseaRes <- fgseaRes[padj < 0.05, ]
sig_fgseaRes <- sig_fgseaRes[order(pval), ]
sig_fgseaRes

# create bar chart of normalised enrichment scores (NES) for top10 hits
# hidden to avoid overwriting
# png(file = 'plots/FGSEA_KANSL1_mut_vs_WT_excl_BaSq_padj0.05.png', width = 12, height = 5, units = 'in', res = 1000)

ggplot(sig_fgseaRes, aes(x = NES, y=reorder(pathway, -pval), fill = factor(sign(NES)))) + 
  geom_bar(stat = "identity", width = 0.8) +
  labs(title = "GSEA", x = "Normalised Enrichment Score (NES)", y = "Pathway") +
  theme_minimal(base_size = 16) +
  scale_fill_manual(values = c("#0754A2", "#B10029"), guide = "none") +
  scale_y_discrete(labels = function(x) gsub("^HALLMARK_", "", x)) +
  theme(axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# dev.off()

# keratinisation geneset still here, as well as a gene set related to metastasis in melanoma
# this could support the idea that mutating KANSL1 has some sort of role in promoting EMT/metastasis

# extract leading edge genes for each of these genesets
krt_ledge <- unlist(sig_fgseaRes[pathway == 'REACTOME_KERATINIZATION', leadingEdge])
met_ledge <- unlist(sig_fgseaRes[pathway == 'WINNEPENNINCKX_MELANOMA_METASTASIS_UP', leadingEdge])
krt_ledge
met_ledge

# save data
gsea_df <- as.data.frame(sig_fgseaRes)
gsea_df$leadingEdge <- sapply(gsea_df$leadingEdge, paste, collapse = ",")
write.csv(gsea_df, "results/FGSEA_KANSL1_muts_vs_WT_excl_BaSq.csv", row.names = FALSE)

# 9.3 Clinical Stage ------------------------------------------------------

# idea: mutations in KANSL1 might facilitate dedifferentiation, EMT and metastasis
# are KANSL1 mutants enriched for stage IV cancers?
# want to look at luminal cancers only (no Ba/Sq) 
# use exclusions from 9.2.1

clinical_metadata <- read.table("data/TCGA-BLCA-clinical-metadata.tsv", 
                            header = TRUE, sep = "\t")
summary(con_class)

# add Patient ID column with 01A removed from each of the IDs
clinical_metadata$Patient_ID <- substr(clinical_metadata$cases.submitter_id, 1, 12)

# what are the pathological stages for the KANSL1 mutants we extracted?
path_stage_kansl1_mut <- clinical_metadata[clinical_metadata$Patient_ID %in% kansl1_muts, ] |> 
  select('Patient_ID', 'diagnoses.ajcc_pathologic_stage')

# add column saying TRUE/FALSE for whether each cancer is stage IV
path_stage_kansl1_mut$StageIV <- FALSE
path_stage_kansl1_mut$StageIV <- grepl("Stage IV", path_stage_kansl1_mut$diagnoses.ajcc_pathologic_stage) # grepl() returns true if a thing is inside another thing
head(path_stage_kansl1_mut)

# what are the pathological stages for the KANSL1 WTs we extracted, after Ba/Sq tumours
path_stage_kansl1_wt <- clinical_metadata[clinical_metadata$Patient_ID %in% kansl1_wt, ] |> 
  select('Patient_ID', 'diagnoses.ajcc_pathologic_stage')

# add column saying TRUE/FALSE for whether each cancer is stage IV
path_stage_kansl1_wt$StageIV <- FALSE
path_stage_kansl1_wt$StageIV <- grepl("Stage IV", path_stage_kansl1_wt$diagnoses.ajcc_pathologic_stage) # grepl() returns true if a thing is inside another thing
head(path_stage_kansl1_wt)

# look at stages using pie charts

# set up plot matrix
par(mfrow = c(1,2))

# pie chart of entire cohort
StageIV_counts_wt <- table(path_stage_kansl1_wt$StageIV)
pie(StageIV_counts_wt,
    labels = names(StageIV_counts_wt),
    main = 'Stage IV counts of entire BLCA TCGA cohort excl Ba/Sq')

# pie chart of KANSL1 mutants
StageIV_counts_kansl1_mut <- table(path_stage_kansl1_mut$StageIV)
pie(StageIV_counts_kansl1_mut,
    labels = names(StageIV_counts_kansl1_mut),
    main = 'Stage IV counts of KANSL1 mutants')

# no, KANSL1 mutants are not enriched for stage IV tumours
# observationally, there is the same proportion of Stage IV tumours