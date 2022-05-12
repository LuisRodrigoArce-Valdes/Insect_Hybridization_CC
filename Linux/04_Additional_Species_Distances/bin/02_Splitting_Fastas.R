#!/usr/bin/Rscript
# 02_Splitting_Fastas.R
# Made by Luis Rodrigo Arce Valdés, to filter fasta files, choose COI sequences only and split per species
# Calling up libraries
library(phylotools)
library(dplyr)
library(tidyr)
library(stringr)
rm(list = ls())

# Groups vector
groups <- c("Extra")

# Creating list
fastas <- list()

# List of expected species per group
species <- list()
for (i in groups) {
  species[[i]] <- read.table(paste0("../results/",i,"/01_BOLD/species.txt"), sep = ",", header = F, col.names = "Species")
}

# Reading fasta files files per order
for (i in groups) {
  fastas[[i]] <- read.fasta(paste0("../results/",i,"/01_BOLD/seqs.fasta"), clean_name = F)
  # Filtering for COI
  fastas[[i]] <- fastas[[i]][grepl("COI-5P", fastas[[i]]$seq.name),]
  # Splitting column names
  fastas[[i]] <- separate(data = fastas[[i]], col = seq.name, into = c("BOLD", "Species", "Gene", "GenBank"), sep = "\\|", fill = "right")
  # Now we will remove "subspecies", to do that we will split the first appearance of " " by a "_":
  fastas[[i]]$Species <- sub(" ","_",fastas[[i]]$Species)
  # Now removing all text after a space, i. e. subspecies epitet
  fastas[[i]]$Species <- sub(" .*", "", fastas[[i]]$Species)
  # Sorting by species names
  fastas[[i]] <- fastas[[i]][order(fastas[[i]]$Species, decreasing = F),]
  row.names(fastas[[i]]) <- 1:nrow(fastas[[i]])
  # Saving table (to reference in suplementary material)
  write.table(fastas[[i]], paste0("../results/",i,"/01_BOLD/info.tsv"), sep = "\t", quote = F, row.names = F)
  # Looking for missing species
  missings <- setdiff(species[[i]]$Species, sub("_"," ",fastas[[i]]$Species))
  # Writing missing species
  write(missings, paste0("../results/",i,"/01_BOLD/missings.txt"), ncolumns = 1)
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
