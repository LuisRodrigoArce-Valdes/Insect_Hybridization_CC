#!/usr/bin/Rscript
# 02_Corrected_Genetic_Distancess.R
# Made by Luis Rodrigo Arce Valdés, to estimate genetic distances between COIs of hybridising species, using phylogenetic correction
rm(list = ls())

# Calling up libraries
library(ape)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggtext)

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

# Saving table
consensus <- consensus[order(consensus$Order, consensus$Cross),]
hybrids <- consensus
hybrids <- separate(hybrids, col=Cross, into=c("Sp1","Sp2"), sep = "_X_")
hybrids$Sp1 <- gsub("_"," ", hybrids$Sp1)
hybrids$Sp2 <- gsub("_"," ", hybrids$Sp2)
write.table(hybrids, "../figures/02_Pairs_Distances.tsv", quote = F, sep = "\t", row.names = F)
rm(hybrids)

# Phylogenetic correction####
# To phylogenetically correct we need to average comparsions with multiple species

# Adding again to a list (for for-looping)
distances <- list()
distances[["consensus"]] <- consensus
#distances[["longest"]] <- longest
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
#longest <- distances$longest
rm(distances, means, i, order, s, sp2, u, sp.order)

# Plots and Statistics ####
# New directory
dir.create("../figures/", showWarnings = F)

# Reading tables with meta information
# List of hybridizing species per group
species <- list()
for (i in groups) {
  species[[i]] <- read.delim(paste0("../../01_Consensus/data/",i,"_Hybrids.txt"), header = T, fileEncoding = "ISO-8859-1")
  species[[i]]$Order <- i
}

# Merging dataframes
species <- bind_rows(species)
species <- species[,c(10,1:9)]

# Formatting for supplementary material
# First ordering alphabetically
species <- species[order(species$Order, species$Family, species$Sp1, species$Sp2),]

# Creating reference column
species$Authors <- gsub("\\.$","",species$Authors)
species$Authors <- gsub(" *$","",species$Authors)
species$Authors <- gsub("\\(","",species$Authors)
species$Authors <- gsub("\\)","",species$Authors)
species$Reference <- paste0(species$Authors,", ",species$Year)
species$Reference <- gsub(", NA","",species$Reference)

# Printing table (easier to read than multiple excels)
write.table(species[-c(2,7,8,9,10)], "../figures/01_Hybrids_Cases_Table.tsv", quote = F, row.names = F, sep = "\t")

# Creating cross column
species$Cross <- paste0(species$Sp1, "_X_", species$Sp2)
species$Cross <- gsub(" ","_", species$Cross)

# Now we will keep the "consensus" strategy and add meta information columns
# COIs <- COIs[COIs$Strategy=="Consensus",-4]
COIs <- consensus
rm(consensus)
COIs$Sp1 <- gsub("_"," ", COIs$Sp1)
COIs$Sp2 <- gsub("_"," ", COIs$Sp2)

# And adding meta info
Suborder <- vector()
Family <- vector()

# Widening info
species <- gather(species, "Position", "Species", 4:5)

# Removing white spaces at the end of species names
species$Species <- gsub(" *$","",species$Species)

for (i in 1:nrow(COIs)) {
  sp <- COIs[i,"Sp1"]
    Suborder <- append(Suborder, species[species$Species==sp,"Suborder"][1])
    Family <- append(Family, species[species$Species==sp,"Family"][1])
}

# Tiny fix
Condition[Condition=="FIeld"] <- "Field"

# Adding columns
COIs <- data.frame(COIs, Suborder=Suborder, Family=Family)

# Ordering columns
COIs <- COIs[,c(1,5,6,2,3,4)]
COIs <- COIs[order(COIs$Order, COIs$Family, COIs$Sp1, COIs$Sp2),]

# Writing
write.table(COIs, "../figures/03_Corrected_Pairs_Distances.tsv", quote = F, sep = "\t", row.names = F)

# Plotting
colors <- c('#66c2a5','#fc8d62','#8da0cb','#e78ac3','#a6d854')
colors2 <- c('#7fc97f','#beaed4','#fdc086','#ffff99','#386cb0','#f0027f','#bf5b17','#666666')

# Violin plot per order
png("../figures/04_Violins.png", width = 2400, height = 1600)
ggplot(COIs) +
  geom_violin(aes(y=Distance, x=Order, fill=Order), scale = "width", draw_quantiles = 0.50) +
  geom_point(aes(y=Distance, x=Order), alpha=0.75, size=6) +
  theme_classic() +
  scale_fill_manual(values = colors) +
  labs(y = "Genetic Distance", caption = "Outlier = *Gryllus texensis* X *Gryllus rubens*")+
  theme(text = element_text(size = 45, family = "serif"),
        legend.position = "none",
        axis.title.x = element_blank(),
        plot.caption = element_markdown())
dev.off()

# Removing outlier
COIs <- COIs[COIs$Distance < 0.30,]

# Plotting again
png("../figures/05_Violins.png", width = 2400, height = 1600)
ggplot(COIs) +
  geom_violin(aes(y=Distance, x=Order, fill=Order), alpha=0.5, scale = "width", draw_quantiles = 0.50) +
  geom_point(aes(y=Distance, x=Order), alpha=0.5, size=6) +
  theme_classic() +
  scale_fill_manual(values = colors) +
  labs(y = "Genetic Distance") +
  theme(text = element_text(size = 45, family = "serif"),
        legend.position = "none",
        axis.title.x = element_blank())
dev.off()

# Statistical testing:
sink("../figures/06_Distances_testing.txt", append = F, split = T)
kruskal.test(Distance ~ Order, data = COIs)
for(i in c("none","bonferroni","BY")) {
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print(pairwise.wilcox.test(COIs$Distance, COIs$Order, p.adjust.method = i))
}
sink()

# Factoring suborder
COIs$Suborder <- factor(COIs$Suborder, levels = c("Brachycera", "Nematocera",
                                                  "\"Ant\"", "\"Other\"",
                                                  "\"Moth\"", "\"Butterfly\"",
                                                  "Caelifera", "Ensifera",
                                                  "Zygoptera"))

# Plotting for suborder
png("../figures/07_Violins_Suborders.png", width = 2400, height = 1600)
ggplot(COIs[COIs$Order!="Odonata",]) +
  facet_wrap(~ Order, ncol = 2, scales = "free_x") +
  geom_violin(aes(y=Distance, x=Suborder, fill=Suborder), alpha=0.5, scale = "width", draw_quantiles = 0.50) +
  geom_point(aes(y=Distance, x=Suborder), alpha=0.5, size=6) +
  scale_fill_manual(values = colors2) +
  theme_classic() +
  labs(y = "Genetic Distance") +
  theme(text = element_text(size = 45, family = "serif"),
        legend.position = "none",
        axis.title.x = element_blank())
dev.off()

# Statistical testing:
sink("../figures/08_Distances_testing.txt", append = F, split = T)
for(n in unique(COIs[COIs$Order!="Odonata","Order"])) {
  print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
  print(n)
  print(wilcox.test(Distance ~ Suborder, data = COIs[COIs$Order==n,], exact=F))
}
sink()

# Plotting for families
png("../figures/09_Violins_Families.png", width = 2400, height = 1600)
ggplot(COIs[COIs$Order!="Odonata",]) +
  facet_wrap(~ Order, nrow = 4, scales = "free") +
  geom_violin(aes(y=Distance, x=Family, fill=Suborder), alpha=0.5, scale = "width", draw_quantiles = 0.50) +
  geom_point(aes(y=Distance, x=Family), size=4, alpha=0.5) +
  scale_fill_manual(values = colors2) +
  theme_classic() +
  labs(y = "Genetic Distance")+
  theme(text = element_text(size = 45, family = "serif"),
        legend.position = "bottom",
        axis.text.x = element_text(size=28, angle = 10, hjust = 1),
        axis.title.x = element_blank())
dev.off()
