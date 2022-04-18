#!/usr/bin/Rscript
# 04_Paired_Fastas.R
# Made by Luis Rodrigo Arce Valdés, to merge into a single fasta file hybridizing species COIs
# Calling up libraries
rm(list = ls())
library(phylotools)
library(tidyr)
library(dplyr)

# Groups vector
groups <- c("Orthoptera", "Lepidoptera", "Diptera", "Hymenoptera")

# List of hybridizing species per group
species <- list()
for (i in groups) {
  species[[i]] <- read.delim(paste0("../data/",i,"_Hybrids.txt"), header = T)
  # Keeping only the species columns
  species[[i]] <- species[[i]][,c("Sp1","Sp2")]
  # Replacing " " by "_"
  species[[i]] <- as.data.frame(apply(species[[i]], 2, function(x) gsub(" ", "_", x)))
  # Adding number of pair
  species[[i]]$pair <- 1:nrow(species[[i]])
  # Tidying
  species[[i]] <- gather(species[[i]], "Number", "Species", 1:2)
}

# Reading fasta files files per order
fastas <- list()
for (i in groups) {
  # Listing all .fasta files per group
  temp = list.files(path = paste0("../results/",i,"/04_Consensus/"), pattern = "*.fasta")
  # Editing names
  temp <- gsub(".fasta","",temp)
  # Reading all .fasta files and assigning them to their name
  for (n in temp) {
    print(n)
    fastas[[i]][[n]] <- read.fasta(paste0("../results/",i,"/04_Consensus/",n,".fasta"), clean_name = F)
  }
  fastas[[i]] <- bind_rows(fastas[[i]], .id = "Species")
  fastas[[i]] <- fastas[[i]][,-2]
}

# Now we will add each species sequence
for (i in groups) {
seqs <- vector()
  for (n in 1:nrow(species[[i]])) {
    sp <- species[[i]][n,"Species"]
    if (sp %in% fastas[[i]]$Species) {
      fa <- fastas[[i]][fastas[[i]]$Species==sp,"seq.text"]
      seqs <- append(seqs, fa)  
    }
    else {
      seqs <- append(seqs, NA)
    }
  }
species[[i]]$COI <- seqs
# Removing species without COIs
species[[i]] <- species[[i]][complete.cases(species[[i]]),]
}

# Writing output fasta file per pair of species
for (i in groups) {
  species[[i]] <- species[[i]][,-2]
  fastas[[i]] <- list()
  dir.create(paste0("../results/",i,"/05_Paired_Fastas"), showWarnings = F)
  # Creating indiviual data.frames per pair of species
  for (n in unique(species[[i]]$pair)) {
    fastas[[i]][[n]] <- species[[i]][species[[i]]$pair==n,-1]
    colnames(fastas[[i]][[n]]) <- c("seq.name","seq.text")
    if (nrow(fastas[[i]][[n]])==2) {
      print(nrow(fastas[[i]][[n]]))
      phylotools::dat2fasta(fastas[[i]][[n]], outfile = paste0("../results/",i,"/05_Paired_Fastas/",fastas[[i]][[n]]$seq.name[1], "_X_", fastas[[i]][[n]]$seq.name[2],".fasta"))
    }
  }
}
