#!/bin/sh
# 03_Alignments_and_Consensus.sh
# By Luis Rodrigo Arce Valdés (28/03/22)
# With this script we will align all COI sequences for each species and then create a Consensus sequence for each one
# For looping across all 4 orders

for i in Hymenoptera Orthoptera Diptera #Lepidoptera Odonata
do
	mkdir -p ../data/$i/03_Muscle
	mkdir -p ../data/$i/04_Consensus
	for n in $(ls ../data/$i/02_Species_Fastas/)
	do
		echo "########"
		echo "$i - $n"
		grep ">" ../data/$i/02_Species_Fastas/$n | wc -l
		# Aligning with muscle's "super5" algorithm (faster than the original method in large databases)
		muscle -super5 ../data/$i/02_Species_Fastas/$n -output ../data/$i/03_Muscle/$n.afa
		# Creating consesus sequence [cons only works if we have more than two sequencues]
		if [ $(grep ">" ../data/$i/02_Species_Fastas/$n | wc -l) -gt 1 ]
		then
			cons -sequence ../data/$i/03_Muscle/$n.afa -outseq ../data/$i/04_Consensus/$n -name $n
		else
			cp ../data/$i/03_Muscle/$n.afa ../data/$i/04_Consensus/$n
			sp=">$n"
			sed -i "1s/.*/$sp/" ../data/$i/04_Consensus/$n
		fi
	done
done
