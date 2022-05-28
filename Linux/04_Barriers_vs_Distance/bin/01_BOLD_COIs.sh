#!/bin/sh
# 01_Getting_COIs.sh
# By Luis Rodrigo Arce Valdés (28/03/22)
# Using this script we will download COI sequences from the BOLD database hybridizing pair species
# http://www.boldsystems.org/
# https://github.com/CNuge/BOLD-CLI

# 01.- Copying odondata species
cp ../../../Windows/Busquedas/Barriers/Odonata/results/odonates.tsv ../data/01_raw/Odonata_Barriers.tsv

# 02.- Copying 4 original orders input files
for i in Lepidoptera Orthoptera Hymenoptera Diptera
do
cp ../../../Windows/Busquedas/Barriers/$i/${i}_Barriers.tsv ../data/01_raw/
done

# For looping across all 5 orders
for i in Odonata Lepidoptera Orthoptera Hymenoptera Diptera
do
echo $i

# 03. Creating output diretories
mkdir -p ../data/$i
mkdir -p ../data/$i/01_BOLD

# 04.- First I am doing a list of all hybridizing species
# First column (First species):
cut -f1 ../data/01_raw/${i}_Barriers.tsv | tail -n +2 > tmp1

# Second column (Second species):
cut -f2 ../data/01_raw/${i}_Barriers.tsv | tail -n +2 > tmp2

# Concatenating both columns, sorting and removing duplicates
cat tmp1 tmp2 | sort | uniq | sed 's/ *$//' > ../data/$i/01_BOLD/${i}_species.txt
rm tmp1 tmp2

# 05.- Searching COI sequences within the BOLD systems database
./bold-cli -output ../data/$i/01_BOLD/$i.fasta -query sequence -marker COI-5P -taxon ../data/$i/01_BOLD/${i}_species.txt

done