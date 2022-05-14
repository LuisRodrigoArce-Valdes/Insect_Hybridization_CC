#!/bin/sh
# 05_Paired_Alignments.sh
# By Luis Rodrigo Arce Valdés (28/03/22)
# With this script we will align all consensus COI sequences per pair of species
# For looping across all 4 orders

for i in Odonata Orthoptera Hymenoptera Diptera #Lepidoptera
do
	mkdir -p ../data/$i/06_Paired_Muscle
	for n in $(ls ../data/$i/05_Paired_Fastas/)
	do
		echo "########"
		echo "$i - $n"
		# Aligning with muscle
		muscle -align ../data/$i/05_Paired_Fastas/$n -output ../data/$i/06_Paired_Muscle/$n
	done
done
