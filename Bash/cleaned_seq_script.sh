#!/bin/bash
# Script for removing taxa without dna sequences for topological test w/ iqtree
# this script must be used on raws data, not aligned sequences

mkdir -p mol_data_cleaned  # create folder if it doesn't exist

for file in *.fasta; do  # remove all ligne with NNN[...]NN AND the ligne above (taxa name)
    sed '$!N;/NNNNNNNNNNNNNNNNNNNNN/!P;D' "$file" > "mol_data_cleaned/${file%.fasta}_cleaned.fasta"
done

echo "Task is completed"