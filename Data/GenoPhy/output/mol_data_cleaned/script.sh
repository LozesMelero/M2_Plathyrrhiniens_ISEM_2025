#!/bin/bash

# Create the new directory if it doesn't exist
mkdir -p alignement_cleaned_data

# Loop through the cleaned files and run MAFFT
for file in *_cleaned.fasta; do
    # Run MAFFT on the file and save the output in the new directory
	mafft --auto "$file" > "alignement_cleaned_data/$(basename "$file" .fasta)_aligned.fasta"
done

echo "Alignment completed."
