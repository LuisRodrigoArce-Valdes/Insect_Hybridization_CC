#!/bin/sh
# 01_Getting_COIs.sh
# By Luis Rodrigo Arce Valdés (28/03/22)
# Using this script we will download COI sequences from the BOLD database hybridizing pair species
# Adapted for additional species
# http://www.boldsystems.org/
# https://github.com/CNuge/BOLD-CLI

# 01. Creating output diretories
mkdir -p ../results/Extra/01_BOLD

# 02.- First I am doing a list of all hybridizing species
# Second column (First species):
cut -f2 ../data/Hybrids.txt > tmp1

# Third column (Second species):
cut -f3 ../data/Hybrids.txt > tmp2

# Concatenating both columns, sorting and removing duplicates
cat tmp1 tmp2 | sort | uniq > ../results/Extra/01_BOLD/species.txt
rm tmp1 tmp2

# 03.- Searching COI sequences within the BOLD systems database
./bold-cli -output ../results/Extra/01_BOLD/seqs.fasta -query sequence -marker COI-5P -taxon ../results/Extra/01_BOLD/species.txt
