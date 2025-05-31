#!/bin/bash
###########################
# Name: local_blast_Mico.sh
# Author: Lucas Buffan
# Aim: Performs local blast search of our nuclear markers in the published genomes of Costa-Araújo et al. (2019)
##########################

# Requires ncbi-blast+
# sudo apt install ncbi-blast+

# Create local BLAST database
#makeblastdb -in ./Mico_database/Mico_db.fasta -dbtype nucl -out ./Mico_database/Mico_db

# BLAST each sequence against the database
for seq in ABCA1 RAG1 ADORA3 ERC2 FES FOXP1 MAPKAP1 RPGRIP1 SIM1 ZFX
do
	blastn -task blastn -query ./Out_queries/$seq/*.fasta -db ./Mico_database/Mico_db -out ./Out_queries/$seq/results_blast_$seq.txt -outfmt 7
done
