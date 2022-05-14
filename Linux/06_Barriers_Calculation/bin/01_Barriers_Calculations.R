#!/usr/bin/Rscript
# 01_Barriers_Calculations.R
# Made by Luis Rodrigo Arce Valdés, to estimate RI using excel files and correlate it with genetic distance
rm(list = ls())

# Calling libraries
library(dplyr)
library(ggplot2)

# Four original orders ####

# Groups vector
groups <- c("Orthoptera", "Diptera", "Hymenoptera") #"Lepidoptera"

# Barriers per group
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
  barriers[[i]][barriers[[i]]$RI=="YYYY","pos"] <- 5
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

# Subseting odonates COIs
COI.odonata <- COI[COI$Order=="Odonata",]

# Merging isolation barriers and genetic distance
COI <- COI[match(intersect(COI$Cross, isolation$Cross),COI$Cross),]
isolation <- isolation[match(intersect(COI$Cross, isolation$Cross), isolation$Cross),]

# Checking if vectors are identical
identical(isolation$Cross, COI$Cross)

# Merging dataframes
isolation <- cbind(isolation, Gen=COI$Distance)
rm(COI)

# Plotting
dir.create("../figures", showWarnings = F)
png("../figures/01_Scatterplot.png", width = 1800, height = 1000, units = "px", res = 300)
ggplot(isolation) +
  geom_point(aes(x=Gen, y=RI, color=Order)) +
  theme_classic()
dev.off()

# Now transforming into categorical variables
# First we will unify for hymenopterans
isolation[isolation==0.25] <- 0.20
isolation[isolation==0.50] <- 0.40
isolation[isolation==0.75] <- 0.60

# Now transforming ri into a factor
isolation <- isolation[order(isolation$RI, decreasing = T),]
isolation$RI <- factor(isolation$RI, levels = unique(isolation$RI),
                       labels = c("Choice_Mating", "Oviposition", "Hybrid_Survival",
                                  "Hybrid_Fertility","Haldanes_Rule","Complete_Gene_Flow"))
# Plotting
png("../figures/02_Barriers.png", width = 2400, height = 1200, units = "px", res = 300)
ggplot(isolation) +
  facet_wrap(. ~ Order, nrow = 3) +
  geom_boxplot(aes(x=RI, y=Gen, color=Order)) +
  geom_point(aes(x=RI, y=Gen, color=Order)) +
  theme_classic()
dev.off()

# Odonates ####

# Now we will include the odonate species
# Reading table of isolation in odonates
odonata <- read.delim("../../05_Barriers_vs_Distance/data/01_raw/Odonata_Barriers.tsv")

# Replacing "-1" by "0"
odonata[odonata==-1] <- 0

# Estimating RI
odonata$RI <- 1-(odonata$RI/6)

# Editing isolation species names
odonata$Cross <- paste0(odonata$Male,"_X_",odonata$Female)
odonata$Cross <- gsub(" ","_", odonata$Cross)
odonata.RI <- odonata[,c(10,9)]

# Merging odonata barriers and genetic distance
COI.odonata <- COI.odonata[match(intersect(COI.odonata$Cross, odonata.RI$Cross),COI.odonata$Cross),]
odonata.RI <- odonata.RI[match(intersect(COI.odonata$Cross, odonata.RI$Cross), odonata.RI$Cross),]

# Checking if vectors are identical
identical(odonata.RI$Cross, COI.odonata$Cross)

# Merging dataframes
odonata.RI <- cbind(odonata.RI, Gen=COI.odonata$Distance)

# Plotting
png("../figures/03_Odonata_Scatterplot.png", width = 1800, height = 1000, units = "px", res = 300)
ggplot(odonata.RI) +
  geom_point(aes(x=Gen, y=RI)) +
  theme_classic()
dev.off()

# Now transforming ri into a factor
odonata.RI <- odonata.RI[order(odonata.RI$RI, decreasing = T),]
odonata.RI$RI <- factor(odonata.RI$RI, levels = unique(odonata.RI$RI),
                       labels = c("Sexual", "TandemAttempt", "Tandem",
                                      "Mating","Oviposition","Hybrid","CompleteGeneFlow"))

# Plotting
png("../figures/04_Odonata_Barriers.png", width = 2400, height = 1200, units = "px", res = 300)
ggplot(odonata.RI) +
  geom_boxplot(aes(x=RI, y=Gen)) +
  geom_point(aes(x=RI, y=Gen)) +
  geom_point(data = odonata.RI[odonata.RI$Cross=="Ischnura_elegans_X_Ischnura_graellsii" | odonata.RI$Cross=="Ischnura_graellsii_X_Ischnura_elegans",], aes(x=RI, y=Gen, color=Cross)) +
  labs(caption = "Se retendrá la cruza recíproca que más anvance en la reproducción. Por ejemplo, Macho elegans con Hembra graellsii") +
  theme_classic() +
  theme(legend.position = "bottom")
dev.off()


# Removing reciprocal crosses in Odonata ####

# To remove reciprocal crosses we will sort cross column alphabetically
odonata$code <- vapply(odonata$Cross, function(xi) paste(sort(strsplit(xi, NULL)[[1]]), collapse=''), '')

# Sorting first by code then by isolation
odonata <- odonata[order(odonata$code, odonata$RI),]

# Removing duplicated codes (species pairs)
odonata <- odonata[!duplicated(odonata$code),]

# Repeating analyses
odonata.RI <- odonata[,c(10,9)]

# Merging odonata barriers and genetic distance
COI.odonata <- COI.odonata[match(intersect(COI.odonata$Cross, odonata.RI$Cross),COI.odonata$Cross),]
odonata.RI <- odonata.RI[match(intersect(COI.odonata$Cross, odonata.RI$Cross), odonata.RI$Cross),]

# Checking if vectors are identical
identical(odonata.RI$Cross, COI.odonata$Cross)

# Merging dataframes
odonata.RI <- cbind(odonata.RI, Gen=COI.odonata$Distance)

# Plotting
png("../figures/05_Odonata_Scatterplot_1D.png", width = 1800, height = 1000, units = "px", res = 300)
ggplot(odonata.RI) +
  geom_point(aes(x=Gen, y=RI)) +
  theme_classic()
dev.off()

# Now transforming ri into a factor
odonata.RI <- odonata.RI[order(odonata.RI$RI, decreasing = T),]
odonata.RI$RI <- factor(odonata.RI$RI, levels = unique(odonata.RI$RI),
                        labels = c("TandemAttempt", "Tandem",
                                       "Mating","Oviposition","Hybrid","CompleteGeneFlow"))

# Plotting
png("../figures/06_Odonata_Barriers_1D.png", width = 2400, height = 1200, units = "px", res = 300)
ggplot(odonata.RI) +
  geom_boxplot(aes(x=RI, y=Gen)) +
  geom_point(aes(x=RI, y=Gen)) +
  geom_point(data = odonata.RI[odonata.RI$Cross=="Ischnura_elegans_X_Ischnura_graellsii",], aes(x=RI, y=Gen), color="red") +
  theme_classic() +
  theme(legend.position = "bottom")
dev.off()


# Merging Odonates and other insects ####
odonata.RI$Order <- "Odonata"
odonata.RI <- odonata.RI[,c(4,1,2,3)]

# Binding
isolation <- rbind(isolation, odonata.RI)

# Now unifing names of barriers
barriers <- data.frame(Original=unique(isolation$RI), New=c("Mating","Oviposition","Hybrid",
                                                "Complete_Gene_Flow","Complete_Gene_Flow","Complete_Gene_Flow",
                                                "Mating","Mating","Mating",
                                                "Hybrid","Complete_Gene_Flow"))

# Viewing barriers
barriers

# Factoring
isolation$RI <- factor(isolation$RI, levels = barriers$Original, labels = barriers$New)

# Plotting
# Plotting
png("../figures/07_Barriers_Merged.png", width = 2400, height = 1200, units = "px", res = 300)
ggplot(isolation) +
  facet_wrap(. ~ Order, nrow = 2) +
  geom_boxplot(aes(x=RI, y=Gen, color=Order)) +
  geom_point(aes(x=RI, y=Gen, color=Order)) +
  theme_classic() +
  theme(legend.position = "none")
dev.off()

# Now scaterploting
isolation$RI <- factor(isolation$RI, levels = unique(isolation$RI), labels = c("1","0.666","0.333","0"))
isolation$RI <- as.numeric(as.vector(isolation$RI))

# Plotting
png("../figures/08_Scatterplot_Merged.png", width = 1800, height = 1000, units = "px", res = 300)
ggplot(isolation) +
  geom_point(aes(x=Gen, y=RI, color=Order)) +
  theme_classic()
dev.off()
