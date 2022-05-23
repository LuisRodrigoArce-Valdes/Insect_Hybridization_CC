#!/usr/bin/Rscript
# 02_Corrected_Genetic_Distancess.R
# Made by Luis Rodrigo Arce Valdés, to estimate genetic distances between COIs of hybridising species, using phylogenetic correction
rm(list = ls())

# Calling up libraries
library(ape)
library(ggplot2)
library(dplyr)
library(tidyr)

# Groups vector
groups <- c("Odonata", "Orthoptera", "Lepidoptera", "Diptera", "Hymenoptera")

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
rm(n, temp)

# Merging dataframes
longest <- bind_rows(longest)

# Phylogenetic correction####
# To phylogenetically correct we need to average comparsions with multiple species

# Adding again to a list (for for-looping)
distances <- list()
distances[["consensus"]] <- consensus
distances[["longest"]] <- longest
rm(consensus, longest)

# For looping in both datasets
for (u in 1:length(distances)) {
  # First we will remove repeated comparisons
  distances[[u]] <- separate(distances[[u]], col=Cross, into=c("Sp1","Sp2"), sep = "_X_")
  
  sp.order <- data.frame()
  # Sorting species
  for (i in 1:nrow(distances[[u]])) {
    sp.order[i,c(1,2)] <- sort(as.character(distances[[u]][i,c(2,3)]))
  }
  
  # Adding columns to distances[[u]]
  distances[[u]][,c(2,3)] <- sp.order[,c(1,2)]
  
  # Removing duplicates
  distances[[u]] <- unique(distances[[u]])
  
  # Ordering per species
  distances[[u]] <- distances[[u]][order(distances[[u]]$Sp1, distances[[u]]$Sp2),]
  
  # Averaging repeated species
  means <- vector()
  sp2 <- vector()
  order <- vector()
  
  # For looping
  for (i in unique(distances[[u]]$Sp1)) {
    p <- mean(as.numeric(distances[[u]][distances[[u]]$Sp1==i,4]))
    if (length(distances[[u]][distances[[u]]$Sp1==i,4])==1) {
      s <- distances[[u]][distances[[u]]$Sp1==i,3]
    }
    else {
      s <- paste0("AAA",i,"_average")
    }
    order <- append(order, distances[[u]][distances[[u]]$Sp1==i,1][1])
    means <- append(means, p)
    sp2 <- append(sp2, s)
  }
  
  # Dataframing
  distances[[u]] <- data.frame(Order=order,
                          Sp1=unique(distances[[u]]$Sp1),
                          Sp2=sp2,
                          Distance=means)
  
  # Sorting by second species
  distances[[u]] <- distances[[u]][order(distances[[u]]$Sp2, distances[[u]]$Sp1),c(1,3,2,4)]
  
  # And repeating averages
  sp.order <- data.frame()
  for (i in 1:nrow(distances[[u]])) {
    sp.order[i,c(1,2)] <- sort(as.character(distances[[u]][i,c(2,3)]))
  }
  
  # Adding columns to distances[[u]]
  distances[[u]][,c(3,2)] <- sp.order[,c(1,2)]
  
  # Ordering per species
  distances[[u]] <- distances[[u]][order(distances[[u]]$Sp2, distances[[u]]$Sp1),]
  
  # Averaging repeated species
  means <- vector()
  sp2 <- vector()
  order <- vector()
  
  # For looping
  for (i in unique(distances[[u]]$Sp2)) {
    p <- mean(as.numeric(distances[[u]][distances[[u]]$Sp2==i,4]))
    if (length(distances[[u]][distances[[u]]$Sp2==i,4])==1) {
      s <- distances[[u]][distances[[u]]$Sp2==i,3]
    }
    else {
      s <- "Average"
    }
    order <- append(order, distances[[u]][distances[[u]]$Sp2==i,1][1])
    means <- append(means, p)
    sp2 <- append(sp2, s)
  }
  
  # Dataframing
  distances[[u]] <- data.frame(Order=order,
                          Sp1=unique(distances[[u]]$Sp2),
                          Sp2=sp2,
                          Distance=means)
  
  # Gsubing
  distances[[u]]$Sp2 <- sub("^AAA.*","Average", distances[[u]]$Sp2)
}

# Extracting
consensus <- distances$consensus
longest <- distances$longest
rm(distances, means, i, order, s, sp2, u, sp.order)

# Plots and Statistics ####
consensus$Strategy <- "Consensus"
longest$Strategy <- "Longest"
COIs <- rbind(consensus, longest)
rm(consensus, longest)

# Pasting species names
COIs$Cross <- paste0(COIs$Sp1, "_X_",COIs$Sp2)
COIs <- COIs[,c(1,6,4,5)]

# Widening
Wide.COIs <- spread(COIs, Strategy, Distance)

# Looking for outliers
Wide.COIs$Outliers <- Wide.COIs$Consensus > 0.20 | Wide.COIs$Longest > 0.20

# New directory
dir.create("../figures_corrected/", showWarnings = F)

# Scatterplot
png("../figures_corrected/01_Strategies.png", width = 1800, height = 1600, res = 300)
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
png("../figures_corrected/02_Violin_Outliers.png", width = 1200, height = 800)
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

# With the consensus strategy there seems to not be outliers, we wont delete them since we will be using this strat

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
write.table(species, "../figures_corrected/Merged_Table.tsv", quote = F, row.names = F, sep = "\t")

# Creating cross column
species$Cross <- paste0(species$Sp1, "_X_", species$Sp2)
species$Cross <- gsub(" ","_", species$Cross)

# Now we will keep the "consensus" strategy and add meta information columns
COIs <- COIs[COIs$Strategy=="Consensus",-4]

# To look for the meta info of each species we will split the column again
COIs <- separate(COIs,Cross,c("Sp1","Sp2"),"_X_")
COIs$Sp1 <- gsub("_"," ",COIs$Sp1)

# And adding meta info
Condition <- vector()
Suborder <- vector()
Family <- vector()

# Widening info
species <- gather(species, "Position", "Species", 4:5)

for (i in 1:nrow(COIs)) {
  sp <- COIs[i,"Sp1"]
    Condition <- append(Condition, species[species$Species==sp,"Condition"][1])
    Suborder <- append(Suborder, species[species$Species==sp,"Suborder"][1])
    Family <- append(Family, species[species$Species==sp,"Family"][1])
}

# Tiny fix
Condition[Condition=="FIeld"] <- "Field"

# Adding columns
COIs <- data.frame(COIs, Suborder=Suborder, Family=Family, Condition=Condition)

# Returning to Cross
COIs$Sp1 <- gsub(" ","_", COIs$Sp1)
COIs$Cross <- paste0(COIs$Sp1, "_X_", COIs$Sp2)

# Ordering columns
COIs <- COIs[,c(1,5,6,8,4)]
COIs <- COIs[order(COIs$Order, COIs$Suborder, COIs$Family, COIs$Cross),]

# Also exporting table
write.table(COIs, "../figures_corrected/Distances.txt", quote = F, sep = "\t", row.names = F)

# Filtering odonates
odonates <- COIs[COIs$Order=="Odonata",]
COIs <- COIs[COIs$Order!="Odonata",]

# Factoring suborder
COIs$Suborder <- factor(COIs$Suborder, levels = c("Brachycera", "Nematocera",
                                                  "\"Ant\"", "\"Other\"",
                                                  "\"Moth\"", "\"Butterfly\"",
                                                  "Caelifera", "Ensifera"))

# Plotting for suborder
png("../figures_corrected/03_Violin_Suborder.png", width = 1200, height = 800)
ggplot(COIs) +
  facet_wrap(~ Order, ncol = 2, scales = "free_x") +
  geom_violin(aes(y=Distance, x=Suborder, fill=Suborder), alpha=0.5, scale = "width") +
  geom_jitter(aes(y=Distance, x=Suborder), size=2.5) +
  scale_shape_manual(values = c(1, 3)) +
  theme_classic() +
  labs(y = "Genetic Distance", caption = paste0("Consensus Strategy;  Genetic distance: ", model))+
  theme(text = element_text(size = 20, family = "serif"),
        legend.position = "bottom",
        axis.title.x = element_blank())
dev.off()

# Plotting for families
png("../figures_corrected/04_Violin_Families.png", width = 1200, height = 800)
ggplot(COIs[COIs$Order!="Lepidoptera",]) +
  facet_wrap(~ Order, nrow = 3, scales = "free") +
  geom_violin(aes(y=Distance, x=Family, fill=Suborder), alpha=0.5, scale = "width") +
  geom_jitter(aes(y=Distance, x=Family), size=2.5) +
  scale_shape_manual(values = c(1, 3)) +
  theme_classic() +
  labs(y = "Genetic Distance", caption = paste0("Consensus Strategy;  Genetic distance: ", model))+
  theme(text = element_text(size = 20, family = "serif"),
        legend.position = "bottom",
        axis.title.x = element_blank())
dev.off()

# Butterflies
png("../figures_corrected/05_Violin_Butterflies.png", width = 1200, height = 800)
ggplot(COIs[COIs$Order=="Lepidoptera",]) +
  facet_wrap(~ Suborder, nrow = 2, scales = "free") +
  geom_violin(aes(y=Distance, x=Family, fill=Suborder), alpha=0.5, scale = "width") +
  geom_jitter(aes(y=Distance, x=Family), size=2.5) +
  scale_shape_manual(values = c(1, 3)) +
  theme_classic() +
  labs(y = "Genetic Distance", caption = paste0("Consensus Strategy;  Genetic distance: ", model))+
  theme(text = element_text(size = 20, family = "serif"),
        legend.position = "bottom",
        axis.title.x = element_blank())
dev.off()

