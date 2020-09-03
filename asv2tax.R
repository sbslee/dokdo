library(dada2)

args <- commandArgs(trailingOnly=TRUE)
asv_file <- args[[1]]
taxa_file <- args[[2]]
species_file <- args[[3]]
output_file <- args[[4]]

seqtab.nochim <- as.matrix(read.csv(asv_file, row.names=1))

taxa <- assignTaxonomy(seqtab.nochim, taxa_file, multithread=TRUE)
taxa <- addSpecies(taxa, species_file)

write.csv(taxa, output_file)