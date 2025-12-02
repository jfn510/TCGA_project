# Capstone project analysis

# Set-up ------------------------------------------------------------------

# load libraries
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

# log into google drive
googledrive::drive_auth()
# click 1

# Identify tumours with KANSL1 mutation -----------------------------------

# read Whole Exome Sequencing maf
wxs_maf <- read.maf('data/WXS.maf')

# get summaries of all mutated genes
gene_tots <- getGeneSummary(wxs_maf)
gene_tots[grep('KANSL1', gene_tots$Hugo_Symbol),]$MutatedSamples
gene_tots[grep('KANSL1', gene_tots$Hugo_Symbol),]$AlteredSamples
# not sure what the differnce between the mutated and altered samples columns are
# 24 patients had mutations in KANSL1

# lollipop plot could not be produced, cBioPortal will be used instead

# extract KANSL1 mutant patient IDs
kansl1_muts <- unique(wxs_maf@data[Hugo_Symbol == "KANSL1", Patient_Id])
kansl1_muts


# Investigating consensus classifiers of KANSL1 mutants -------------------


# are the KANSL1 mutants more frequently in any one of the consensus classifiers?
con_class <- read.table("data/mRNA_gc47-TPMs_ConsensusClassifier.tsv", 
                        header = TRUE, sep = "\t")
summary(con_class)

# add Patient ID column with 01A removed from each of the IDs
con_class$Patient_ID <- substr(con_class$ID, 1, 12)

# what are the classifiers for the KANSL1 mutants we extracted?
con_class_kansl1 <- con_class[con_class$Patient_ID %in% kansl1_muts, ] |> 
  select('Patient_ID', 'consensusClass')

# Create pie chart --------------------------------------------------------

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

# DE analysis - KANSL1 mutants vs ALL non-mutants -----------------------------

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

# load RNAseq count data, adapt column names (remove last four characters)
counts <- read.table('data/mRNA_gc47-counts.tsv', check.names = FALSE,
                      header = TRUE, row.names = 1, sep = '\t')

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

# DE analysis - 24vs24 ----------------------------------------------------

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
# 138

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

# Savepoint A -------------------------------------------------------------

# create save point so you can come back to the same random 24
# save.image("savepointA_DEA1-complete.RData")

load("savepointA_DEA1-complete.RData")

# DEA excluding some genes ------------------------------------------------

# create data frame which VEP and SIFT scores can be added to
mut_info <- unique(wxs_maf@data[Hugo_Symbol == "KANSL1", Patient_Id, Protein_Change])
mut_info <- mut_info[order(mut_info$Patient_Id), ] # order Patient IDs alphabetically
mut_info$VEP <- 