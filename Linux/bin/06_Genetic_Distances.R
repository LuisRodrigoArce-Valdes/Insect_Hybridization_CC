#!/usr/bin/Rscript
# 06_Genetic_Distancess.R
# Made by Luis Rodrigo Arce Valdés, to estimate genetic distances between COIs of hybridising species
rm(list = ls())

# Calling up libraries
library(ape)
library(ggplot2)

# Groups vector
groups <- c("Orthoptera", "Lepidoptera", "Diptera", "Hymenoptera")

# Reading FASTAs
COIs <- list()
DIST <- list()

for (i in groups) {
  COIs[[i]] <- list()
  DIST[[i]] <- list()
  # Listing all .fasta files per group
  temp <- list.files(path = paste0("../results/",i,"/06_Paired_Muscle/"), pattern = "*.fasta")
  # Editing names
  temp <- gsub(".fasta","",temp)
  # Reading all .fasta files and assigning them to their name
  for (n in temp) {
    print(n)
    COIs[[i]][[n]] <- read.FASTA(paste0("../results/",i,"/06_Paired_Muscle/",n,".fasta"), type = "DNA")
    DIST[[i]][[n]] <- dist.dna(COIs[[i]][[n]], model = "TN93", variance = F, as.matrix = T)
    DIST[[i]][[n]] <- DIST[[i]][[n]][1,2]
  }
  DIST[[i]] <- as.data.frame(t.data.frame(as.data.frame(DIST[[i]])))
  DIST[[i]]$Cross <- row.names(DIST[[i]])
  row.names(DIST[[i]]) <- 1:nrow(DIST[[i]])
  DIST[[i]] <- DIST[[i]][,c(2,1)]
  colnames(DIST[[i]])[2] <- "Distance"
  DIST[[i]]$Order <- i
  DIST[[i]] <- DIST[[i]][,c(3,1,2)]
}

# Plotting
ggplot(DIST$Orthoptera) +
  geom_violin(aes(y=Distance, x=Order)) +
  geom_jitter(aes(y=Distance, x=Order))
