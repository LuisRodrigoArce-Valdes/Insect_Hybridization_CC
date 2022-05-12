#!/usr/bin/Rscript
# 01_Excel_To_Values.R
# Made by Luis Rodrigo Arce Valdés, to estimate RI using excel files and correlate it with genetic distance
rm(list = ls())

# Calling libraries
library(dplyr)
library(ggplot2)

# Groups vector
groups <- c("Orthoptera", "Diptera", "Hymenoptera") #"Odonata", "Lepidoptera"

# Barrieres per group
barriers <- list()
for (i in groups) {
  barriers[[i]] <- read.table(paste0("../../05_Barriers_vs_Distance/data/01_raw/",i,"_Barriers.tsv"), sep = "\t", header = T)
}

# While estimating reproductive barriers we will do a different estimation for Hymenoptera (because they create hybrids of only one sex)
groups <- groups[groups!="Hymenoptera"]

# Estimating reproductive barriers
for (i in groups) {
  # Adding a column pasting the strings of each row
  barriers[[i]]$RI <- paste0(barriers[[i]]$Choice.Mating,
                             barriers[[i]]$Hybrid.Egg,
                             barriers[[i]]$Hybrid.Adult,
                             barriers[[i]]$One.Sex.Fertile.Hybrids,
                             barriers[[i]]$Fertile.Hybrids)
  # Looking for the last barrier each species can do gene flow
  RI <- vector()
  for (n in 1:nrow(barriers[[i]])) {
    pos <- unlist(gregexpr('N', barriers[[i]]$RI[n]))[1]
    RI <- append(RI, pos)
  }
  barriers[[i]]$pos <- RI
  # Removing -1s that have "?"
  barriers[[i]][barriers[[i]]$RI=="YYYYY","pos"] <- 6
  barriers[[i]] <- barriers[[i]][barriers[[i]]$pos!=-1,]
  # Estimating reproductive isolation barriers
  barriers[[i]]$RI <- 1-((barriers[[i]]$pos-1)/5)
  barriers[[i]] <- barriers[[i]][,-10]
}


# Now for hymenopterans with a small change (RI divided by 4 total barriers)
for (i in "Hymenoptera") {
  # Adding a column pasting the strings of each row
  barriers[[i]]$RI <- paste0(barriers[[i]]$Choice.Mating,
                             barriers[[i]]$Hybrid.Egg,
                             barriers[[i]]$Hybrid.Adult,
                             barriers[[i]]$One.Sex.Fertile.Hybrids)
  # Looking for the last barrier each species can do gene flow
  RI <- vector()
  for (n in 1:nrow(barriers[[i]])) {
    pos <- unlist(gregexpr('N', barriers[[i]]$RI[n]))[1]
    RI <- append(RI, pos)
  }
  barriers[[i]]$pos <- RI
  # Removing -1s that have "?"
  barriers[[i]][barriers[[i]]$RI=="YYYYNA","pos"] <- 5
  barriers[[i]] <- barriers[[i]][barriers[[i]]$pos!=-1,]
  # Estimating reproductive isolation barriers
  barriers[[i]]$RI <- 1-((barriers[[i]]$pos-1)/4)
  barriers[[i]] <- barriers[[i]][,-10]
}

# Creating groups vector
groups <- c(groups, "Hymenoptera")

# Adding column name
for (i in groups) {
  barriers[[i]]$Order <- i
}

# Merging dataframes
isolation <- bind_rows(barriers)

# Ordering columns
isolation <- isolation[,c(10,1,2,9)]

# Reading genetic distance data frame
COI <- read.delim("../../05_Barriers_vs_Distance/results/01_Genetic_Distances.tsv")

# Editing isolation species names
isolation$Cross <- paste0(isolation$Sp1,"_X_",isolation$Sp2)
isolation$Cross <- gsub(" ","_", isolation$Cross)
isolation <- isolation[,c(1,5,4)]

# Merging isolation barriers and genetic distance
COI <- COI[match(intersect(COI$Cross, isolation$Cross),COI$Cross),]
isolation <- isolation[match(intersect(COI$Cross, isolation$Cross), isolation$Cross),]

# Checking if vectors are identical
identical(isolation$Cross, COI$Cross)

# Merging dataframes
isolation <- cbind(isolation, Gen=COI$Distance)
rm(COI)

# Plotting
ggplot(isolation) +
  geom_point(aes(x=Gen, y=RI, color=Order)) +
  theme_classic()
