# Made by Luis Rodrigo Arce Valdés, to filter fasta files and choose COI sequencues only

# Calling up libraries
library(phylotools)
library(dplyr)
library(tidyr)
library(stringr)
rm(list = ls())

# Groups vector
groups <- c("Orthoptera")

# Creating list
fastas <- list()

# List of expected species per group
species <- list()
for (i in groups) {
  species[[i]] <- read.table(paste0("../results/",i,"/01_BOLD/",i,"_species.txt"), sep = ",", header = F, col.names = "Species")
}

# Reading fasta files files per order
for (i in groups) {
  fastas[[i]] <- read.fasta(paste0("../results/",i,"/01_BOLD/",i,".fasta"), clean_name = F)
  # Filtering for COI
  fastas[[i]] <- fastas[[i]][grepl("COI-5P", fastas[[i]]$seq.name),]
  # Splitting column names
  fastas[[i]] <- separate(data = fastas[[i]], col = seq.name, into = c("BOLD", "Species", "Gene", "GenBank"), sep = "\\|", fill = "right")
  # Counting number of non "gapped" characters per fasta sequence
  fastas[[i]]$Nucleotides <- nchar(fastas[[i]]$seq.text) - str_count(fastas[[i]]$seq.text, "-")
  # Sorting by number of nucleotides
  fastas[[i]] <- fastas[[i]][order(fastas[[i]]$Nucleotides, decreasing = T),]
  # Filtering for the first appearance of each species (longest sequenced COI)
  fastas[[i]] <- fastas[[i]][match(unique(fastas[[i]]$Species), fastas[[i]]$Species),]
  # Sorting by species names
  fastas[[i]] <- fastas[[i]][order(fastas[[i]]$Species, decreasing = F),]
  row.names(fastas[[i]]) <- 1:nrow(fastas[[i]])
  # Saving table
  write.table(fastas[[i]], paste0("../results/",i,"/01_BOLD/FinalCOIs_",i,".tsv"), sep = "\t", quote = F, row.names = F)
  # Looking for missing species
  missings <- setdiff(species[[i]]$Species, fastas[[i]]$Species)
  # Writing missing species
  write(missings, paste0("../results/",i,"/01_BOLD/Missings_",i,".txt"), ncolumns = 1)
  # Transforming table to phylotols data
  fastas[[i]] <- data.frame(seq.name=paste(fastas[[i]]$BOLD, fastas[[i]]$Species, fastas[[i]]$Gene, fastas[[i]]$GenBank, sep = "|"), seq.text=fastas[[i]]$seq.text)
  # Writting output fasta
  phylotools::dat2fasta(fastas[[i]], outfile = paste0("../results/",i,"/01_BOLD/FinalCOIs_",i,".fasta"))
}
