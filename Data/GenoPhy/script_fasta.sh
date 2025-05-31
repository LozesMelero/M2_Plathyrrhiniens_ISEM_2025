#!/bin/bash

mkdir -p genophy_cleaned


exclure_taxons=() # Initialiser le tableau vide pour exclure les taxons


while IFS= read -r taxon; do # Lire les taxons dans liste.txt et les ajouter au tableau
    exclure_taxons+=("$taxon")
done < liste.txt


for file in output/*.fasta; do
    cp "$file" "genophy_cleaned/$(basename "$file")"

   
    for taxon in "${exclure_taxons[@]}"; do  # Boucle pour exclure chaque taxon de la liste
        awk -v taxon=">$taxon" '
        BEGIN { in_taxon = 0 }
        {
            if ($0 == taxon) {
                in_taxon = 1
                next
            }
            if (in_taxon && substr($0,1,1) == ">") {
                in_taxon = 0
            }
            if (!in_taxon) {
                print $0
            }
        }' "genophy_cleaned/$(basename "$file")" > temp_file && mv temp_file "genophy_cleaned/$(basename "$file")"
    done

done

echo "Task is completed"