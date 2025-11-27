# to create list of HGV of KANSL1 mutations, to go into Ensembl's VEP
# to then look at SIFT and PolyPhen scores

library(data.table)

# the tsv of KANSL1 mutations was downloaded from cBioPortal
mutations <- fread('data/KANSL1_mutations.tsv')

# assign HGV values to vector
HGVSg <- mutations$HGVSg
HGVSc <- mutations$HGVSc

# write txt
writeLines(HGVSg, 'outputs/HGVSg.txt')
writeLines(HGVSc, 'outputs/HGVSc.txt')
