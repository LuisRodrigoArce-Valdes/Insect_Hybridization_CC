#!/usr/bin/Rscript
# 06_Genetic_Distancess.R
# Made by Luis Rodrigo Arce Valdés, to estimate genetic distances between COIs of hybridising species
rm(list = ls())

# Calling up libraries
library(ape)
library(ggplot2)
library(dplyr)
library(tidyr)

# Groups vector
groups <- c("Extra")

# Genetic Distance Model (help: dist.dna):
model <- "raw"

# Consensus ####
# Reading FASTAs
COIs <- list()
consensus <- list()

for (i in groups) {
  COIs[[i]] <- list()
  consensus[[i]] <- list()
  # Listing all .fasta files per group
  temp <- list.files(path = paste0("../results/",i,"/06_Paired_Muscle/"), pattern = "*.fasta")
  # Editing names
  temp <- gsub(".fasta","",temp)
  # Reading all .fasta files and assigning them to their name
  for (n in temp) {
    print(n)
    COIs[[i]][[n]] <- read.FASTA(paste0("../results/",i,"/06_Paired_Muscle/",n,".fasta"), type = "DNA")
    consensus[[i]][[n]] <- dist.dna(COIs[[i]][[n]], model = model, variance = F, as.matrix = T)
    consensus[[i]][[n]] <- consensus[[i]][[n]][1,2]
  }
  consensus[[i]] <- as.data.frame(t.data.frame(as.data.frame(consensus[[i]])))
  consensus[[i]]$Cross <- row.names(consensus[[i]])
  row.names(consensus[[i]]) <- 1:nrow(consensus[[i]])
  consensus[[i]] <- consensus[[i]][,c(2,1)]
  colnames(consensus[[i]])[2] <- "Distance"
  consensus[[i]]$Order <- i
  consensus[[i]] <- consensus[[i]][,c(3,1,2)]
}

# Merging dataframes
consensus <- bind_rows(consensus)

# Writing results
write.table(consensus, "../results/Extra/Distances.tsv", sep = "\t", quote = F, row.names = F)
