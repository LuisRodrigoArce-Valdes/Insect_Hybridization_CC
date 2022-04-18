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
groups <- c("Orthoptera", "Lepidoptera", "Diptera", "Hymenoptera")

# Genetic Distance Model (help: dist.dna):
model <- "TN93"

# Consensus ####
# Reading FASTAs
COIs <- list()
consensus <- list()

for (i in groups) {
  COIs[[i]] <- list()
  consensus[[i]] <- list()
  # Listing all .fasta files per group
  temp <- list.files(path = paste0("../../01_Consensus/results/",i,"/06_Paired_Muscle/"), pattern = "*.fasta")
  # Editing names
  temp <- gsub(".fasta","",temp)
  # Reading all .fasta files and assigning them to their name
  for (n in temp) {
    print(n)
    COIs[[i]][[n]] <- read.FASTA(paste0("../../01_Consensus/results/",i,"/06_Paired_Muscle/",n,".fasta"), type = "DNA")
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

# Longest ####
# Reading FASTAs
COIs <- list()
longest <- list()

for (i in groups) {
  COIs[[i]] <- list()
  longest[[i]] <- list()
  # Listing all .fasta files per group
  temp <- list.files(path = paste0("../../02_Longest/results/",i,"/04_Paired_Muscle/"), pattern = "*.fasta")
  # Editing names
  temp <- gsub(".fasta","",temp)
  # Reading all .fasta files and assigning them to their name
  for (n in temp) {
    print(n)
    COIs[[i]][[n]] <- read.FASTA(paste0("../../02_Longest/results/",i,"/04_Paired_Muscle/",n,".fasta"), type = "DNA")
    longest[[i]][[n]] <- dist.dna(COIs[[i]][[n]], model = model, variance = F, as.matrix = T)
    longest[[i]][[n]] <- longest[[i]][[n]][1,2]
  }
  longest[[i]] <- as.data.frame(t.data.frame(as.data.frame(longest[[i]])))
  longest[[i]]$Cross <- row.names(longest[[i]])
  row.names(longest[[i]]) <- 1:nrow(longest[[i]])
  longest[[i]] <- longest[[i]][,c(2,1)]
  colnames(longest[[i]])[2] <- "Distance"
  longest[[i]]$Order <- i
  longest[[i]] <- longest[[i]][,c(3,1,2)]
}

# Merging dataframes
longest <- bind_rows(longest)

# Plots and Statistics ####
consensus$Strategy <- "Consensus"
longest$Strategy <- "Longest"
COIs <- rbind(consensus, longest)
rm(consensus, longest)

# Widening
Wide.COIs <- spread(COIs, Strategy, Distance)

# Looking for outliers
Wide.COIs$Outliers <- Wide.COIs$Consensus > 0.25 | Wide.COIs$Longest > 0.25

# Scatterplot
png("../figures/01_Strategies.png", width = 1600, height = 1600, res = 300)
ggplot(Wide.COIs) +
  scale_y_continuous(limits = c(0,1)) +
  scale_x_continuous(limits = c(0,1)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_point(aes(x=Consensus, y=Longest, color=Order)) +
  geom_label(data = Wide.COIs[Wide.COIs$Outliers==T,], 
             aes(x=Consensus, y=Longest, label=Cross),
             nudge_x = 0.1, nudge_y = 0.03, size=2) +
  labs(caption = paste0("Genetic distance: ", model)) +
  theme_classic() +
  theme(text = element_text(size = 12, family = "serif"),
        legend.position = "bottom")
dev.off()

# Violin plots
png("../figures/02_Violin_Outliers.png", width = 1200, height = 800)
ggplot(COIs) +
  facet_wrap(~ Strategy, ncol = 2) +
  geom_violin(aes(y=Distance, x=Order, fill=Order), alpha=0.5) +
  geom_jitter(aes(y=Distance, x=Order), size=2.5) +
  theme_classic() +
  labs(y = "Genetic Distance", caption = paste0("Genetic distance: ", model))+
  theme(text = element_text(size = 20, family = "serif"),
        legend.position = "none",
        axis.title.x = element_blank())
dev.off()

# Without outliers
Not.Outliers <- Wide.COIs[Wide.COIs$Outliers==F,"Cross"]
COIs <- COIs[COIs$Cross %in% Not.Outliers,]

for (i in unique(COIs$Strategy)) {
  png(paste0("../figures/03_",i,"_Without_Outliers.png"), width = 1200, height = 800)
  plot <- ggplot(COIs[COIs$Strategy==i,]) +
    geom_violin(aes(y=Distance, x=Order, fill=Order), alpha=0.5) +
    geom_jitter(aes(y=Distance, x=Order), size=2.5) +
    theme_classic() +
    labs(subtitle = paste0("Strategy: ", i) ,y = "Genetic Distance", caption = paste0("Outliers removed; genetic distance: ", model)) +
    theme(text = element_text(size = 20, family = "serif"),
          legend.position = "none",
          axis.title.x = element_blank())
  print(plot)
  dev.off()
}
