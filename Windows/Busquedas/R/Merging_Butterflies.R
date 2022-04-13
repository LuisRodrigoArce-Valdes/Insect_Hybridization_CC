rm(list = ls())
# Calling up libraries
library(tidyr)

# This script merges my hybrid butterflies database with that of Pregraves, 2002.
butter <- read.csv("Lepidoptera_Hybrids.csv", header = T)
Pres <- read.csv("Tabla_Presgraves_2002.csv", header = T)

# To look for non-overlapping species pairs first we need to sort alphabetically each species at each pair
# Adding species identification row
butter$num <- 1:nrow(butter)
Pres$num <- 1:nrow(Pres)

# Tidying
butter <- gather(butter, "Pair", "Species", 1:2)
Pres <- gather(Pres, "Pair", "Species", 2:3)

# Sorting alphabetically
butter <- butter[order(butter$num, butter$Species),]
Pres <- Pres[order(Pres$num, Pres$Species),]

# Editing species 1 and species 2
butter$Pair <- rep(c("Sp1", "Sp2"), nrow(butter)/2)
Pres$Pair <- rep(c("Sp1", "Sp2"), nrow(Pres)/2)

# Widening
butter <- spread(butter, Pair, Species)
Pres <- spread(Pres, Pair, Species)

# Creating merged species column
butter$merged <- paste0(butter$Sp1,"-",butter$Sp2)
Pres$merged <- paste0(Pres$Sp1,"-",Pres$Sp2)

# Looking for new species pairs in Pres
new <- setdiff(Pres$merged, butter$merged)

# Subsetting Pres for new species
Pres <- Pres[Pres$merged %in% new,]

# Sorting dataframes by species number
butter <- butter[order(butter$num),]
Pres <- Pres[order(Pres$num),]

# Selecting columns to export
butter <- butter[,c(7,8,1:5)]

# Editing Pres condition
Pres[Pres=="Allopatric"] <- "Laboratory"

# Now for Pres species
Pres <- data.frame(Sp1=Pres$Sp1, Sp2=Pres$Sp2, Condition=Pres$Natural, MainTopic="Yes", Authors="Presgraves, DC.", Year=2002, Article="Patterns of postzygotic isolation in Lepidoptera")

# rbinding
butter <- rbind(butter, Pres)

# Exporting
write.table(butter, "butter.tsv", quote = F, row.names = F, sep = "\t")
