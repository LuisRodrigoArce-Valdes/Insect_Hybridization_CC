# Script to noquote out title, authors, years and abstract of all found manuscripts
rm(list = ls())

# Select Web  of Science table to analyse
file <- paste0("../Barriers/Odonata2/Odonata.txt")

# Output file name
output <- gsub(".txt","", file)

# Add search information
engine <- "Web of Science"
string <- '(AB=("reproductive barriers" OR "reproductive isolation") OR AK=("reproductive barriers" OR "reproductive isolation") OR TI=("reproductive barriers" OR "reproductive isolation")) AND (ALL=(Odonata))'
date <- "23/05/2022"


# Reading file (coding for this tab delimited text files is VERY strange)
head <- read.table(file, header = F, nrows = 1, sep = "\t")
file <- read.table(file, header = F, sep = "\t", skip = 1, fill = T, quote = "")

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
