# Data formating
# Using this script we will format Rosas databases to be as the ones we created for RI barriers
rm(list = ls())
library(dplyr)

# Reading files
groups <- c("S1","S2","S3")
odonates <- list()
for (i in groups) {
  odonates[[i]] <- read.delim(paste0("../data/",i,".txt"))
}

# Editing S1 table to be as the other ones
odonates$S1$Genera.1 <- odonates$S1$Genera
odonates$S1 <- odonates$S1[,c(1,2,3,12,4:11)]

# Now editing as groups
for (i in groups) {
  colnames(odonates[[i]]) <- c("Table","Genera.male","Species.male","Genera.female","Species.female",
                                    "Sexual","TandemAttempt","Tandem","Mating","Oviposition","Hybrid.Field","Hybrid.Lab")
}

# It will be easier to merge the tables
odonates <- bind_rows(odonates)

# species epithet to lower
odonates$Species.male <- tolower(odonates$Species.male)
odonates$Species.female <- tolower(odonates$Species.female)

# Removing white spaces that could be at the end of each species
odonates$Genera.male <- gsub(" *$","", odonates$Genera.male)
odonates$Species.male <- gsub(" *$","", odonates$Species.male)
odonates$Genera.female <- gsub(" *$","", odonates$Genera.female)
odonates$Species.female <- gsub(" *$","", odonates$Species.female)

# Removing white spaces that could be at the beggining of each species epithet
odonates$Species.male <- gsub("^ *","", odonates$Species.male)
odonates$Species.female <- gsub("^ *","", odonates$Species.female)

# Removing subspecies epithet
odonates$Species.male <- sub(" .*", "", odonates$Species.male)
odonates$Species.female <- sub(" .*", "", odonates$Species.female)

# Merging male and female columns
odonates$Male <- paste(odonates$Genera.male, odonates$Species.male)
odonates$Female <- paste(odonates$Genera.female, odonates$Species.female)

# Ordering columns
odonates <- odonates[,c(13,14,6:12)]

# Replacing NA
odonates[is.na(odonates)] <- ""

# Pasting Hybrid columns
odonates$Hybrid <- paste0(odonates$Hybrid.Field, odonates$Hybrid.Lab)

# Removing extra columns
odonates <- odonates[,c(-8,-9)]

# Replacing "" by a character
odonates[odonates==""] <- "N"

# Creating string column of RI
odonates$RI <- paste0(odonates$Sexual, odonates$TandemAttempt, odonates$Tandem,
                     odonates$Mating, odonates$Oviposition, odonates$Hybrid)

# For each row identifying on which column reproduction stops
ri <- vector()
for (i in 1:nrow(odonates)) {
  pos <- tail(unlist(gregexpr('x', odonates$RI[i])), n=1)
  print(i)
  print(pos)
  ri <- append(ri, pos)
}
odonates$RI <- ri

# Saving table (later we will unify it with the other orders)
dir.create("../results", showWarnings = F)
write.table(odonates,"../results/odonates.tsv", sep = "\t", quote = F, row.names = F)
