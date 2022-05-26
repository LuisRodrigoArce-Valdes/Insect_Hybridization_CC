#!/usr/bin/Rscript
# 01_Barriers_Calculations.R
# Made by Luis Rodrigo Arce Valdés, to estimate RI using excel files and correlate it with Genetic Distance
rm(list = ls())

# Calling libraries
library(dplyr)
library(ggplot2)
library(tidyr)
library(scales)
library(ggrepel)

# Reading reproductive isolation tables ####
# Groups vector
groups <- c("Orthoptera", "Diptera", "Hymenoptera", "Lepidoptera")

# Barriers per group
barriers <- list()
for (i in groups) {
  barriers[[i]] <- read.table(paste0("../../04_Barriers_vs_Distance/data/01_raw/",i,"_Barriers.tsv"), sep = "\t", header = T)
  barriers[[i]]$Order <- i
  barriers[[i]] <- barriers[[i]][,c(11,1:10)]
}

# Binding rows
orders <- bind_rows(barriers)
rm(barriers)

# Reading table of isolation in odonates
odonata <- read.delim("../../04_Barriers_vs_Distance/data/01_raw/Odonata_Barriers.tsv")

# Replacing "-1" by "0"
odonata[odonata==-1] <- 0

# Estimating RI
odonata$RI <- 1-(odonata$RI/6)

# Editing isolation species names
odonata$Cross <- paste0(odonata$Male,"_X_",odonata$Female)
odonata$Cross <- gsub(" ","_", odonata$Cross)

# Removing reciprocal crosses in Odonata ####
# To remove reciprocal crosses we will sort cross column alphabetically
odonata$code <- vapply(odonata$Cross, function(xi) paste(sort(strsplit(xi, NULL)[[1]]), collapse=''), '')

# Sorting first by code then by isolation
odonata <- odonata[order(odonata$code, odonata$RI),]

# Removing duplicated codes (species pairs)
odonata <- odonata[!duplicated(odonata$code),]

# Couplying odonata dataset into orders ####
# First we will change the names of the barriers
colnames(orders)[5:9] <- c("Assortative_Mating",
                           "Oviposition",
                           "Hybrid_Inviability",
                           "Haldanes_Rule",
                           "Fertile_F1_Hybrids")

# For hymenopterans we will move their measurements, now being Haldanes Rule = NA
orders[orders$Order=="Hymenoptera","Fertile_F1_Hybrids"] <- orders[orders$Order=="Hymenoptera","Haldanes_Rule"]
orders[orders$Order=="Hymenoptera","Haldanes_Rule"] <- NA

# Now unifying in Odonates
odonata <- odonata[,-c(10:11)]

# Factoring odonata RI into orders categories
odonata <- odonata[order(odonata$RI, decreasing = T),]
odonata$RI <- as.character(factor(odonata$RI, levels = unique(odonata$RI), labels = c("Complete RI",
                                                                         "Complete RI",
                                                                         "Complete RI",
                                                                         "Complete RI",
                                                                         "Assortative_Mating",
                                                                         "Oviposition",
                                                                         "Hybrid_Inviability")))

# Odonata equivalency table
odonata.eq <- data.frame(RI=c("Complete RI", "Assortative_Mating", "Oviposition", "Hybrid_Inviability"),
                         labels=c("N_N_N_NA_NA","Y_N_N_NA_NA","Y_Y_N_NA_NA", "Y_Y_Y_NA_NA"))

# Creating new odonata dataframe
bars <- vector()
codes <- vector()
for(i in 1:nrow(odonata)) {
  bar <- odonata[i,"RI"]
  code <- odonata.eq[odonata.eq$RI==bar,"labels"]
  bars <- append(bars, bar)
  codes <- append(codes, code)
}

# Merging at new data.frame
odonata <- data.frame(odonata[,c(1,2)], RI=bars, Labels=codes)
rm(i, bar, bars, code, codes, odonata.eq)

# Splitting labels
odonata <- separate(odonata, col=Labels, into=colnames(orders)[5:9], sep = "_")

# Uniforming data.tables
orders <- orders[,-4]
orders$Authors <- gsub("\\.$",",",orders$Authors)
orders$Reference <- paste(orders$Authors, orders$Year)
orders <- orders[,-c(9,10)]

# Odonates
odonata$Order <- "Odonata"
odonata$Reference <- NA
odonata <- odonata[,-3]
odonata <- odonata[,c(8,1:7,9)]

# Merging data.frames
colnames(odonata) <- colnames(orders)
orders <- rbind(orders, odonata)
rm(odonata)

# Sorting
orders <- orders[order(orders$Order, orders$Sp1, orders$Sp2),]
orders[orders=="NA"] <- NA

# Writing table
dir.create("../figures", showWarnings = F)
write.table(orders, "../figures/01_Merged_Table.tsv", sep = "\t", quote = F, row.names = F)

# Now removing reference column
orders <- orders[,-9]

# Changing NAs to O
orders[is.na(orders)] <- "O"

# Estimating reproductive barriers by adding a column pasting the strings of each row
orders$RI <- paste0(orders$Assortative_Mating,
                    orders$Oviposition,
                    orders$Hybrid_Inviability,
                    orders$Haldanes_Rule,
                    orders$Fertile_F1_Hybrids)

# Splitting into orders:
normal <- orders[orders$Order!="Odonata" & orders$Order!="Hymenoptera",]
odonata <- orders[orders$Order=="Odonata",]
hymenoptera <- orders[orders$Order=="Hymenoptera",]

# "Normal" orders
# Looking for the first barrier that blocks hybridization
RI <- vector()
for (i in 1:nrow(normal)) {
  pos <- unlist(gregexpr('N', normal[i,"RI"]))[1]
  RI <- append(RI, pos)
}
normal$pos <- RI
# Removing -1s that have "?"
normal[normal$RI=="YYYYY","pos"] <- 6
normal <- normal[normal$pos!=-1,]
# Estimating reproductive isolation barriers
normal$RI <- 1-((normal$pos-1)/5)
normal <- normal[,-10]

# Now for Hymenopterans
# Looking for the first barrier that blocks hybridization
RI <- vector()
for (i in 1:nrow(hymenoptera)) {
  pos <- unlist(gregexpr('N', hymenoptera[i,"RI"]))[1]
  RI <- append(RI, pos)
}
hymenoptera$pos <- RI

# Removing -1s that have "?"
hymenoptera[hymenoptera$RI=="YYYOY","pos"] <- 6
hymenoptera <- hymenoptera[hymenoptera$pos!=-1,]

# Estimating reproductive isolation barriers
hymenoptera$RI <- 1-((hymenoptera$pos-1)/5)
hymenoptera <- hymenoptera[,-10]


# Finally for Odonata
# Looking for the first barrier that blocks hybridization
RI <- vector()
for (i in 1:nrow(odonata)) {
  pos <- unlist(gregexpr('N', odonata[i,"RI"]))[1]
  RI <- append(RI, pos)
}
odonata$pos <- RI

# We don't have "?"
odonata[odonata$RI=="YYYOO","pos"] <- 6

# Estimating reproductive isolation barriers
odonata$RI <- 1-((odonata$pos-1)/5)
odonata <- odonata[,-10]

# Merging dataframes
orders <- rbind(normal, hymenoptera, odonata)
rm(normal, hymenoptera, odonata)

# Replacing "O" by NA 
orders[orders=="O"] <- NA

# Adding column for RI Barrier and writing table
orders <- orders[order(orders$RI, decreasing = T),]
orders$Barrier <- factor(orders$RI, levels = unique(orders$RI), labels = c(colnames(orders[4:6]), "Hybrid Sterility", colnames(orders[7:8])))
orders <- orders[order(orders$Order, orders$Sp1, orders$Sp2),]
orders$Barrier <- gsub("_"," ", orders$Barrier)
write.table(orders, "../figures/02_Reproductive_Isolation.tsv", sep = "\t", quote = F, row.names = F)

# Removing useless columns
orders <- orders[,-c(4:8)]

# Bars plots, sorting by RI and factorin barriers
orders <- orders[order(orders$RI, decreasing = F),]
orders[orders=="Haldanes Rule"] <- "Haldane's Rule"
orders$Barrier <- factor(orders$Barrier, levels = unique(orders$Barrier))

# Creating table
isolation.table <- as.data.frame(prop.table(table(orders$Barrier, orders$Order), margin = 2))
colnames(isolation.table)[1:2] <- c("Barrier","Order")

# Widening table for labels
wide.table <- spread(isolation.table, Order, Freq)
wide.table <- wide.table[order(nrow(wide.table):1),]
row.names(wide.table) <- wide.table$Barrier
wide.table <- wide.table[,-1]

# Cumsuming columns
wide.table <- as.data.frame(apply(wide.table, 2, cumsum))

# Returning to tidy table
wide.table$Barrier <- row.names(wide.table)
wide.table <- gather(wide.table, Order, Freq, 1:5)

# Sorting tables
isolation.table$Barrier <- as.character(isolation.table$Barrier)
isolation.table <- isolation.table[order(isolation.table$Barrier, isolation.table$Order),]
wide.table <- wide.table[order(wide.table$Barrier, wide.table$Order),]
isolation.table$CumSum <- wide.table$Freq
rm(wide.table)
isolation.table$Barrier <- factor(isolation.table$Barrier, levels = unique(orders$Barrier))
isolation.table <- isolation.table[order(isolation.table$Barrier),]

# Plotting
png("../figures/03_Bar-Plots.png", width = 24, height = 12, units = "cm", res = 300)
ggplot(isolation.table) +
  geom_col(aes(x=Order, y=Freq, fill=Barrier)) +
  geom_label_repel(data=isolation.table[isolation.table$Freq!=0,],
                   aes(x=Order, y=CumSum, label=paste0(round(Freq*100,1), "%")), 
                   family="serif", size=3, force=0.001, force_pull = 10) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c('#edf8fb','#bfd3e6','#9ebcda','#8c96c6','#8856a7','#810f7c')) +
  labs(y="Frequency") +
  theme_classic() +
  theme(text = element_text(family = "serif"),
        axis.title.x = element_blank())
dev.off()
rm(isolation.table)

# Reading Genetic Distance data frame
COI <- read.delim("../../04_Barriers_vs_Distance/results/01_Genetic_Distances.tsv")

# Editing orders species names
orders$Cross <- paste0(orders$Sp1,"_X_",orders$Sp2)
orders$Cross <- gsub(" ","_", orders$Cross)

# Merging orders barriers and Genetic Distance
COI <- COI[match(intersect(COI$Cross, orders$Cross),COI$Cross),]
orders <- orders[match(intersect(COI$Cross, orders$Cross), orders$Cross),]

# Checking if vectors are identical
identical(orders$Cross, COI$Cross)

# Merging dataframes
orders <- cbind(orders, Gen=COI$Distance)
rm(COI)
orders <- orders[,-6]
orders <- orders[,c(1,2,3,4,6,5)]

# Saving table before phylogenetic correction
write.table(orders, "../figures/04_Iso_Gen.tsv", sep="\t", quote = F, row.names = F)

# I decided to do 2 boxplots, one before and one after phylogenetic correction
# Removing outliers
tmp <- orders[orders$Gen < 0.40,]
tmp$Barrier <- as.character(tmp$Barrier)
tmp <- tmp[order(tmp$RI, decreasing = T),]
tmp$Barrier <- factor(tmp$Barrier, levels = unique(tmp$Barrier))

# Plotting
colors <- c('#66c2a5','#fc8d62','#8da0cb','#e78ac3','#a6d854')

png("../figures/10_Barriers_Non_Corrected.png", width = 12, height = 24, units = "cm", res = 300)
ggplot(tmp) +
  facet_wrap(. ~ Order, ncol = 1, scales = "free_y") +
  geom_boxplot(aes(x=Barrier, y=Gen, fill=Order)) +
  geom_point(aes(x=Barrier, y=Gen), alpha=0.5) +
  scale_fill_manual(values = colors) +
  labs(y="Genetic Distance") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 7, angle = 10, hjust = 1),
        axis.title.x = element_blank(),
        text = element_text(family = "serif"),
        strip.text = element_text(color = "white", size = 0),
        legend.position = "none")
dev.off()

# Statistical testing
sink("../figures/11_Barriers_NonCorrected_Testing.txt", append = F, split = T)
for(i in unique(tmp$Order)) {
  print("###########################")
  print(i)
  print(kruskal.test(Gen ~ Barrier, data = tmp[tmp$Order==i,]))
  for(c in c("none","bonferroni","BY")) {
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print(pairwise.wilcox.test(tmp[tmp$Order==i,"Gen"], tmp[tmp$Order==i,"Barrier"], p.adjust.method = c))
  }
}
sink()
rm(tmp)

# Phylogenetic correction ####
# First we will remove repeated comparisons
sp.order <- data.frame()
# Sorting species
for (i in 1:nrow(orders)) {
  sp.order[i,c(1,2)] <- sort(as.character(orders[i,c(2,3)]))
}

# Adding columns to orders
orders[,c(2,3)] <- sp.order[,c(1,2)]

# Removing duplicates
orders <- unique(orders)

# Ordering per species
orders <- orders[order(orders$Sp1, orders$Sp2),]

# Averaging repeated species
RI <- vector()
Gen <- vector()
sp2 <- vector()
order <- vector()

# For looping
for (i in unique(orders$Sp1)) {
  ri <- mean(as.numeric(orders[orders$Sp1==i,"RI"]))
  gen <- mean(as.numeric(orders[orders$Sp1==i,"Gen"]))
  if (length(orders[orders$Sp1==i,4])==1) {
    s <- orders[orders$Sp1==i,3]
  }
  else {
    s <- paste0("AAA",i,"_average")
  }
  RI <- append(RI, ri)
  order <- append(order, orders[orders$Sp1==i,1][1])
  Gen <- append(Gen, gen)
  sp2 <- append(sp2, s)
}

# Dataframing
orders <- data.frame(Order=order,
                             Sp1=unique(orders$Sp1),
                             Sp2=sp2,
                             RI=RI,
                             Gen=Gen)
rm(gen, Gen, i, order, pos, ri, RI, s, sp2)

# Sorting by second species
orders <- orders[order(orders$Sp2, orders$Sp1),c(1,3,2,4,5)]

# And repeating averages
sp.order <- data.frame()
for (i in 1:nrow(orders)) {
  sp.order[i,c(1,2)] <- sort(as.character(orders[i,c(2,3)]))
}

# Adding columns to orders
orders[,c(3,2)] <- sp.order[,c(1,2)]

# Ordering per species
orders <- orders[order(orders$Sp2, orders$Sp1),]

# Averaging repeated species
RI <- vector()
Gen <- vector()
sp2 <- vector()
order <- vector()

# For looping
for (i in unique(orders$Sp2)) {
  ri <- mean(as.numeric(orders[orders$Sp2==i,"RI"]))
  gen <- mean(as.numeric(orders[orders$Sp2==i,"Gen"]))
  if (length(orders[orders$Sp2==i,4])==1) {
    s <- orders[orders$Sp2==i,3]
  }
  else {
    s <- "Average"
  }
  RI <- append(RI, ri)
  order <- append(order, orders[orders$Sp2==i,1][1])
  Gen <- append(Gen, gen)
  sp2 <- append(sp2, s)
}

# Dataframing
orders <- data.frame(Order=order,
                        Sp1=unique(orders$Sp2),
                        Sp2=sp2,
                        RI=RI,
                        Gen=Gen)
rm(gen, Gen, i, order, ri, RI, s, sp2)

# Gsubing
orders$Sp2 <- sub("^AAA.*","Average", orders$Sp2)

# Scatterplots ####
# Function to get p value from linnear regressions
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

# Global regression
reg <- lm(formula = RI ~ Gen, data = orders)
reg <- paste0("r2 = ", round(summary(reg)$r.squared, 5), "; p = ", round(lmp(reg), 3))

# Plotting
png("../figures/05_Insects.png", width = 18, height = 10, units = "cm", res = 300)
ggplot(orders) +
  geom_smooth(aes(x=Gen, y=RI), formula = y~x, method = "lm", se = F) +
  labs(caption = paste0("Outlier: Enallagma signatum X Enallagma pollutum; ",reg)) +
  geom_point(aes(x=Gen, y=RI)) +
  theme_classic()
dev.off()

# Removing outlier:
# Global regression
orders <- orders[orders$Gen < 0.30,]
reg <- lm(formula = RI ~ Gen, data = orders)
sink("../figures/06_Insects_Correlation.txt", split = F)
print("~~~~~~~~~~~~~~~~~~~~~~~~~~+")
print("ALL INSECTS")
summary(reg)
sink()

# Plotting
png("../figures/07_Insects_Scatterplot.png", width = 18, height = 10, units = "cm", res = 300)
ggplot(orders) +
  geom_point(aes(x=Gen, y=RI, color=Order), alpha=0.75) +
  scale_color_manual(values = colors) +
  labs(x="Genetic Distance", y="Reproductive Isolation") +
  theme_classic() +
  theme(text = element_text(family = "serif"),
        legend.position = "bottom")
dev.off()

# Per orders linear regressions
sink("../figures/08_Orders_Regressions.txt", split = T)
for (i in unique(orders$Order)){
  print("~~~~~~~~~~~~~~~~~~~~~~~~~~~")
  print(i)
  x <- lm(formula = RI ~ Gen, data = orders[orders$Order==i,])
  print(summary(x))
}
sink()

# Plotting
png("../figures/09_Orders_Scatterplots.png", width = 10, height = 20, units = "cm", res = 300)
ggplot(orders) +
  facet_wrap(. ~ Order, ncol = 1) +
  geom_point(aes(x=Gen, y=RI, color=Order)) +
  scale_color_manual(values = colors) +
  labs(x="Genetic Distance", y="Reproductive Isolation") +
  theme_classic() +
  theme(legend.position = "none",
        text = element_text(family = "serif", size=14),
        strip.text = element_text(color = "white", size = 0))
dev.off()

# Now transforming into categorical variables
# We will round RI its closest 0.20
orders$Barrier <- (2*round(orders$RI*10/2))/10

# Now transforming ri into a factor
orders <- orders[order(orders$RI, decreasing = T),]
orders$Barrier <- factor(orders$Barrier, levels = unique(orders$Barrier),
                       labels = c("Assortative Mating", "Oviposition", "Hybrid Inviability",
                                  "Hybrid Sterility","Haldane's Rule","Fertile F1 Hybrids"))

png("../figures/12_Barriers_Corrected.png", width = 12, height = 24, units = "cm", res = 300)
ggplot(orders) +
  facet_wrap(. ~ Order, ncol = 1, scales = "free_y") +
  geom_boxplot(aes(x=Barrier, y=Gen, fill=Order)) +
  geom_point(aes(x=Barrier, y=Gen), alpha=0.5) +
  scale_fill_manual(values = colors) +
  labs(y="Genetic Distance") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 10, angle = 20, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=16),
        text = element_text(family = "serif"),
        strip.text = element_text(color = "white", size = 0),
        legend.position = "none")
dev.off()

# Statistical testing
sink("../figures/13_Barriers_Corrected_Testing.txt", append = F, split = T)
for(i in unique(orders$Order)) {
  print("###########################")
  print(i)
  print(kruskal.test(Gen ~ Barrier, data = orders[orders$Order==i,]))
  for(c in c("none","bonferroni","BY")) {
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print(pairwise.wilcox.test(orders[orders$Order==i,"Gen"], orders[orders$Order==i,"Barrier"], p.adjust.method = c))
  }
}
sink()

