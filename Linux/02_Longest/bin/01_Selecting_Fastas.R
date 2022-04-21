#!/usr/bin/Rscript
# 01_Selecting_Fastas.R
# Made by Luis Rodrigo Arce Valdés, to choose the longest COI sequencue of each species {from the COIs downloaded in the alternative protocol}
# Calling up libraries
library(phylotools)
library(dplyr)
library(tidyr)
library(stringr)
rm(list = ls())

# Groups vector
groups <- c("Odonata", "Orthoptera", "Lepidoptera", "Diptera", "Hymenoptera")

# Creating list
fastas <- list()

# List of expected species per group
species <- list()
for (i in groups) {
  species[[i]] <- read.table(paste0("../../01_Consensus/results/",i,"/01_BOLD/",i,"_species.txt"), sep = ",", header = F, col.names = "Species")
}

# Reading fasta files files per order
for (i in groups) {
  fastas[[i]] <- read.fasta(paste0("../../01_Consensus/results/",i,"/01_BOLD/",i,".fasta"), clean_name = F)
  # Filtering for COI
  fastas[[i]] <- fastas[[i]][grepl("COI-5P", fastas[[i]]$seq.name),]
  # Splitting column names
  fastas[[i]] <- separate(data = fastas[[i]], col = seq.name, into = c("BOLD", "Species", "Gene", "GenBank"), sep = "\\|", fill = "right")
  # Now we will remove "subspecies", to do that we will split the first appearance of " " by a "_":
  fastas[[i]]$Species <- sub(" ","_",fastas[[i]]$Species)
  # Now removing all text after a space, i. e. subspecies epitet
  fastas[[i]]$Species <- sub(" .*", "", fastas[[i]]$Species)
  # Counting number of non "gapped" or "Missing data" characters per fasta sequence
  fastas[[i]]$Nucleotides <- nchar(fastas[[i]]$seq.text) - str_count(fastas[[i]]$seq.text, "-") - str_count(fastas[[i]]$seq.text, "N")
  # Sorting by number of nucleotides
  fastas[[i]] <- fastas[[i]][order(fastas[[i]]$Nucleotides, decreasing = T),]
  # Filtering for the first appearance of each species (longest sequenced COI)
  fastas[[i]] <- fastas[[i]][match(unique(fastas[[i]]$Species), fastas[[i]]$Species),]
  # Sorting by species names
  fastas[[i]] <- fastas[[i]][order(fastas[[i]]$Species, decreasing = F),]
  row.names(fastas[[i]]) <- 1:nrow(fastas[[i]])
  # Saving table (to reference in suplementary material)
  dir.create(paste0("../results/",i,"/01_BOLD"), showWarnings = F, recursive = T)
  write.table(fastas[[i]], paste0("../results/",i,"/01_BOLD/",i,"_info.tsv"), sep = "\t", quote = F, row.names = F)
  # Deleting number of nucleotides column
  fastas[[i]] <- fastas[[i]][,-6]
  # Looking for missing species
  missings <- setdiff(species[[i]]$Species, sub("_"," ",fastas[[i]]$Species))
  # Writing missing species
  write(missings, paste0("../results/",i,"/01_BOLD/",i,"_missings.txt"), ncolumns = 1)
  # Splitting per species
  dir.create(paste0("../results/",i,"/02_Species_Fastas"), showWarnings = F)
  species[[i]] <- list()
  for (n in unique(fastas[[i]]$Species)) {
    species[[i]][[n]] <- fastas[[i]][fastas[[i]]$Species==n,]
    # Transforming table to phylotools data
    species[[i]][[n]] <- data.frame(seq.name=paste(species[[i]][[n]]$BOLD, species[[i]][[n]]$Species, species[[i]][[n]]$Gene, species[[i]][[n]]$GenBank, sep = "|"), seq.text=species[[i]][[n]]$seq.text)
    # Writing output fastas
    phylotools::dat2fasta(species[[i]][[n]], outfile = paste0("../results/",i,"/02_Species_Fastas/",n,".fasta"))
  }
}
