#!/bin/sh
# 03_Alignments_and_Consensus.sh
# By Luis Rodrigo Arce Valdés (28/03/22)
# With this script we will align all COI sequences for each species and then create a Consensus sequence for each one
# For looping across all 4 orders

for i in Hymenoptera Orthoptera Diptera Lepidoptera
do
	mkdir -p ../results/$i/03_Muscle
	mkdir -p ../results/$i/04_Consensus
	for n in $(ls ../results/$i/02_Species_Fastas/)
	do
		echo "########"
		echo "$i - $n"
		grep ">" ../results/$i/02_Species_Fastas/$n | wc -l
		# Aligning with muscle's "super5" algorithm (faster than the original method in large databases)
		muscle -super5 ../results/$i/02_Species_Fastas/$n -output ../results/$i/03_Muscle/$n.afa
		# Creating consesus sequence [cons only works if we have more than two sequencues]
		if [ $(grep ">" ../results/$i/02_Species_Fastas/$n | wc -l) -gt 1 ]
		then
			cons -sequence ../results/$i/03_Muscle/$n.afa -outseq ../results/$i/04_Consensus/$n -name $n
		else
			cp ../results/$i/03_Muscle/$n.afa ../results/$i/04_Consensus/$n
			sp=">$n"
			sed -i "1s/.*/$sp/" ../results/$i/04_Consensus/$n
		fi
	done
done
