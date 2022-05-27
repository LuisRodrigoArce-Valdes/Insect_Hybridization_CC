# 02_COIs_Table.R
# With this script we will create a table with each BOLD and GeneBank record key for each species
rm(list = ls())
library(tidyr)
library(stringr)

# Species included in the document ####
# Reading tables of species used at each analyses
Hybrids <- read.delim("../../02_Genetic_Distances/figures/02_Pairs_Distances.tsv")
Barriers <- read.delim("../../05_Barriers_Calculation/figures/04_Iso_Gen.tsv")

# Tidying tables
Hybrids <- gather(Hybrids[,-4], "Number", "Species", 2:3)
Hybrids <- Hybrids[,-2]

Barriers <- gather(Barriers[,-c(4,5,6)], "Number","Species", 2:3)
Barriers <- Barriers[,-2]

# Merging tables
Species <- rbind(Hybrids, Barriers)
rm(Hybrids, Barriers)

# Removing duplicates
Species <- unique(Species)

# Sorting
Species <- Species[order(Species$Order, Species$Species),]
row.names(Species) <- 1:nrow(Species)

# COIs meta information ####
# BOLD
BOLD <- read.table("../data/Bold.txt", header = F, sep = ">")
BOLD <- data.frame(BOLD[,-1])

# Removing duplicates
BOLD <- unique(BOLD)

# Splitting columns
BOLD <- separate(BOLD,BOLD....1.,c("BOLD","Sp","Gen","GenBank"), sep = "\\|")
BOLD <- BOLD[,c(2,1,4)]
BOLD$Sp <- gsub("_"," ",BOLD$Sp)

# GenBank
GenBank <- read.table("../data/GB.txt", header = F, sep = ">")
GenBank <- data.frame(GenBank[,-1])

# Removing duplicates
GenBank <- unique(GenBank)

# Splitting columns
GenBank <- separate(GenBank,GenBank....1.,as.character(1:5), sep = "\\|")
GenBank <- GenBank[,c(5,4)]
GenBank$`5` <- gsub("^ ","", GenBank$`5`)
GenBank <- separate(GenBank,`5`, c("Gen","sp"), sep = " ", extra = "drop")
GenBank$sp <- paste(GenBank$Gen, GenBank$sp)
GenBank$BOLD <- "NA"
GenBank <- GenBank[,c(2,4,3)]
colnames(GenBank) <- colnames(BOLD)

# Binding data.frames
BOLD <- rbind(BOLD, GenBank)
rm(GenBank)

# Removing species which we have COI but didnt use in the analyses
BOLD <- BOLD[BOLD$Sp %in% Species$Species,]

# For each species at BOLD adding order
orders <- vector()
for(i in 1:nrow(BOLD)) {
  sp <- BOLD[i,"Sp"]
  or <- Species[Species$Species==sp,"Order"]
  orders <- append(orders, or)
}

# Appending
BOLD$Order <- orders

# Ordering
BOLD <- BOLD[,c(4,1,2,3)]
BOLD <- BOLD[order(BOLD$Order, BOLD$Sp),]
row.names(BOLD) <- 1:nrow(BOLD)
colnames(BOLD)[2] <- "Species"
length(unique(BOLD$Species))
BOLD$Sequence <- row.names(BOLD)
BOLD <- BOLD[,c(5,1:4)]

# Exporting
dir.create("../results/", showWarnings = F)
write.table(BOLD, "../results/COIs_Meta_Info.tsv", quote = F, row.names = F, sep = "\t")
