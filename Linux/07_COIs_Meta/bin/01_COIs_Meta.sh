#!/bin/sh
# 01_COIs_Meta.sh
# By Luis Rodrigo Arce Valdés (26/05/22)
# Using this script we will create a database of the meta information of all COI sequences used in this project

# Creating files with the meta information of both BOLD and GenBank sequences
# Hybridization ###
# BOLD
mkdir ../data/
grep "COI-5P" ../../01_Consensus/results/*/02_Species_Fastas/* > ../data/Bold.txt
# GeneBank
grep "|gb|" ../../01_Consensus/results/*/02_Species_Fastas/* > ../data/GB.txt

# Barriers ###
grep "COI-5P" ../../04_Barriers_vs_Distance/data/*a/02_Species_Fastas/* >> ../data/Bold.txt
grep "|gb|" ../../04_Barriers_vs_Distance/data/*a/02_Species_Fastas/* >> ../data/GB.txt



