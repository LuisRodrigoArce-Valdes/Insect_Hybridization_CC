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
model <- "raw"

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
  geom_point(aes(x=Consensus, y=Longest, color=Order), alpha = 0.75) +
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
  geom_violin(aes(y=Distance, x=Order, fill=Order), alpha=0.5, scale = "width") +
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
COIs <- COIs[complete.cases(COIs),]

png("../figures/03_Violin_Not_Outliers.png", width = 1200, height = 800)
ggplot(COIs) +
  facet_wrap(~ Strategy, ncol = 2) +
  geom_violin(aes(y=Distance, x=Order, fill=Order), alpha=0.5, scale = "width") +
  geom_jitter(aes(y=Distance, x=Order), size=2.5) +
  theme_classic() +
  labs(y = "Genetic Distance", caption = paste0("Removed Outliers; Genetic distance: ", model))+
  theme(text = element_text(size = 20, family = "serif"),
        legend.position = "none",
        axis.title.x = element_blank())
dev.off()

# Again the scatterplot
Wide.COIs <- Wide.COIs[Wide.COIs$Outliers==F,]
Wide.COIs <- Wide.COIs[complete.cases(Wide.COIs),]

png("../figures/04_Strategies_Not_Outliers.png", width = 1600, height = 1600, res = 300)
ggplot(Wide.COIs) +
  geom_abline(slope = 1, intercept = 0) +
  geom_point(aes(x=Consensus, y=Longest, color=Order), alpha = 0.75) +
  geom_label(data = Wide.COIs[Wide.COIs$Outliers==T,], 
             aes(x=Consensus, y=Longest, label=Cross),
             nudge_x = 0.1, nudge_y = 0.03, size=2) +
  labs(caption = paste0("Genetic distance: ", model)) +
  theme_classic() +
  theme(text = element_text(size = 12, family = "serif"),
        legend.position = "bottom")
dev.off()

# Reading tables with meta information
# List of hybridizing species per group
species <- list()
for (i in groups) {
  species[[i]] <- read.delim(paste0("../../01_Consensus/data/",i,"_Hybrids.txt"), header = T)
  species[[i]]$Order <- i
}

# Merging dataframes
species <- bind_rows(species)
species <- species[,c(10,1:9)]

# Printing table (easier to read than multiple excels)
write.table(species, "../../Merged_Table.tsv", quote = F, row.names = F, sep = "\t")

# Creating cross column
species$Cross <- paste0(species$Sp1, "_X_", species$Sp2)
species$Cross <- gsub(" ","_", species$Cross)

# Now we will keep the "consensus" strategy and add meta information columns
COIs <- COIs[COIs$Strategy=="Consensus",-4]

# And adding meta info
Condition <- vector()
Suborder <- vector()
Family <- vector()
for (i in 1:nrow(COIs)) {
  sp <- COIs[i,"Cross"]
    Condition <- append(Condition, species[species$Cross==sp,"Condition"][1])
    Suborder <- append(Suborder, species[species$Cross==sp,"Suborder"][1])
    Family <- append(Family, species[species$Cross==sp,"Family"][1])
}

# Tiny fix
Condition[Condition=="FIeld"] <- "Field"

# Adding columns
COIs <- data.frame(COIs, Suborder=Suborder, Family=Family, Condition=Condition)

# Ordering columns
COIs <- COIs[,c(1,4,5,2,6,3)]

# Also exporting table
write.table(COIs, "../../Distances.txt", quote = F, sep = "\t", row.names = F)

# Factoring suborder
COIs$Suborder <- factor(COIs$Suborder, levels = c("Brachycera", "Nematocera",
                                                  "\"Ant\"", "\"Other\"",
                                                  "\"Moth\"", "\"Butterfly\"",
                                                  "Caelifera", "Ensifera"))

# Plotting for Condition
png("../figures/05_Violin_Condition.png", width = 1200, height = 800)
ggplot(COIs) +
  facet_wrap(~ Order, ncol = 2) +
  geom_violin(aes(y=Distance, x=Condition, color=Condition), scale = "width") +
  geom_jitter(aes(y=Distance, x=Condition, fill=Suborder), size=2.5, shape=21) +
  theme_classic() +
  labs(y = "Genetic Distance", caption = paste0("Consensus Strategy; Removed Outliers; Genetic distance: ", model))+
  theme(text = element_text(size = 20, family = "serif"),
        legend.position = "bottom",
        axis.title.x = element_blank())
dev.off()

# Plotting for suborder
png("../figures/06_Violin_Suborder.png", width = 1200, height = 800)
ggplot(COIs) +
  facet_wrap(~ Order, ncol = 2, scales = "free_x") +
  geom_violin(aes(y=Distance, x=Suborder, fill=Suborder), alpha=0.5, scale = "width") +
  geom_jitter(aes(y=Distance, x=Suborder, shape=Condition), size=2.5) +
  scale_shape_manual(values = c(1, 3)) +
  theme_classic() +
  labs(y = "Genetic Distance", caption = paste0("Consensus Strategy; Removed Outliers; Genetic distance: ", model))+
  theme(text = element_text(size = 20, family = "serif"),
        legend.position = "bottom",
        axis.title.x = element_blank())
dev.off()

# Plotting for families
png("../figures/07_Violin_Families.png", width = 1200, height = 800)
ggplot(COIs[COIs$Order!="Lepidoptera",]) +
  facet_wrap(~ Order, nrow = 3, scales = "free") +
  geom_violin(aes(y=Distance, x=Family, fill=Suborder), alpha=0.5, scale = "width") +
  geom_jitter(aes(y=Distance, x=Family, shape=Condition), size=2.5) +
  scale_shape_manual(values = c(1, 3)) +
  theme_classic() +
  labs(y = "Genetic Distance", caption = paste0("Consensus Strategy; Removed Outliers; Genetic distance: ", model))+
  theme(text = element_text(size = 20, family = "serif"),
        legend.position = "bottom",
        axis.title.x = element_blank())
dev.off()

# Butterflies
png("../figures/08_Violin_Butterflies.png", width = 1200, height = 800)
ggplot(COIs[COIs$Order=="Lepidoptera",]) +
  facet_wrap(~ Suborder, nrow = 2, scales = "free") +
  geom_violin(aes(y=Distance, x=Family, fill=Suborder), alpha=0.5, scale = "width") +
  geom_jitter(aes(y=Distance, x=Family, shape=Condition), size=2.5) +
  scale_shape_manual(values = c(1, 3)) +
  theme_classic() +
  labs(y = "Genetic Distance", caption = paste0("Consensus Strategy; Removed Outliers; Genetic distance: ", model))+
  theme(text = element_text(size = 20, family = "serif"),
        legend.position = "bottom",
        axis.title.x = element_blank())
dev.off()
