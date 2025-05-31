#!/bin/bash

# Create the new directory if it doesn't exist
mkdir -p ../Data/Alignment/Fasta_files

# Loop through the files and run MAFFT
for file in ../Data/Raw/*.fasta; do
    # Run MAFFT on the file and save the output in the new directory
    mafft --auto "$file" > "../Data/Alignment/Fasta_files/$(basename "$file" .fasta)_aligned.fasta"
done

echo "Alignment completed."
