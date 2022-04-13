# Script to noquote out title, authors, years and abstract of all found manuscripts
rm(list = ls())

###############################################################################################################
# EDIT THIS SECTION

# Select Web  of Science table to analyse
file <- "../Hybrid/Hymenoptera/09_Hybrid_Hymenoptera.txt"

# Output file name
output <- gsub(".txt","", file)

# Add search information
engine <- "Web of Science"
string <- '(AB=(hybrid OR hybridization OR hybridisation) OR AK=(hybrid OR hybridization OR hybridisation) OR TI=(hybrid OR hybridization OR hybridisation)) AND (AB=(hymenoptera) OR TI=(hymenoptera) or AK=(hymenoptera)) NOT ALL=(fluorescent OR fluorescence OR "in situ" OR resistance)'
date <- "28/03/2022"

###############################################################################################################

# Reading file (coding for this tab delimited text files is VERY strange)
head <- read.table(file, header = F, nrows = 1, sep = "\t")
file <- read.table(file, header = F, sep = "\t", skip = 1, fill = T, quote = "")

# Articles records have an extra tab at the end
file <- file[,-ncol(file)]

# Editing header
head <- as.character(head[1,])
head[1] <- "PT"

# Using as colnames
colnames(file) <- head
rm(head)

# Sorting papers by publication year (we need the first hybridisation records)
# file <- file[order(file$PY, decreasing = F),]

# Adding bold to hybrid
file$AB <- gsub("hybrid","**hybrid**", file$AB)
file$AB <- gsub("Hybrid","**Hybrid**", file$AB)


# printing abstracts
sink(paste0(output, ".md"))

# Additional information
cat(noquote(paste0("###", engine, "\n\n")))
cat(noquote(paste0("*", string, "*\n\n")))
cat(noquote(paste0(date, "\n\n")))

# Printing abstracts
for (i in 1:nrow(file)) {
  cat(noquote(paste0("**", file[i,"TI"],"**\n\n")))
  cat(noquote(paste0(file[i, "AU"], ". ")))
  cat(noquote(paste0(file[i, "PY"], "\n\n")))
  cat(noquote(paste0(file[i, "AB"], "\n\n")))
  cat(noquote("---\n\n"))
}
sink()

# Convert to pdf using https://www.markdowntopdf.com/