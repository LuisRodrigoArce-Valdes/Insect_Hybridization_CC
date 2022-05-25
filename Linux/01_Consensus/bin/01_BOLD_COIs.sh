#!/bin/sh
# 01_Getting_COIs.sh
# By Luis Rodrigo Arce Valdés (28/03/22)
# Using this script we will download COI sequences from the BOLD database hybridizing pair species
# http://www.boldsystems.org/
# https://github.com/CNuge/BOLD-CLI

# For looping across all 4 orders and Odonates

for i in Odonata Orthoptera Lepidoptera Hymenoptera Diptera
do
echo $i

# 01. Creating output diretories
mkdir -p ../results/$i
mkdir -p ../results/$i/01_BOLD

# 02.- First I am doing a list of all hybridizing species
# Third column (First species):
cut -f3 ../data/${i}_Hybrids.txt | tail -n +2 > tmp1

# Fourth column (Second species):
cut -f4 ../data/${i}_Hybrids.txt | tail -n +2 > tmp2

# Concatenating both columns, sorting, removing duplicates and deleting white spaces at the end of each line
cat tmp1 tmp2 | sort | uniq | sed 's/ *$//' > ../results/$i/01_BOLD/${i}_species.txt
rm tmp1 tmp2

# 03.- Searching COI sequences within the BOLD systems database
./bold-cli -output ../results/$i/01_BOLD/$i.fasta -query sequence -marker COI-5P -taxon ../results/$i/01_BOLD/${i}_species.txt

done
