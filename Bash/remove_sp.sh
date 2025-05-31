#!/bin/bash

# Répertoire du script
DIR="$(cd "$(dirname "$0")" && pwd)"

# Chemins relatifs à partir du dossier du script
liste="$DIR/../Data/Liste_sp.txt"
fichier="$DIR/../Analysis/Phylogeny/MrBayes/TED/Matrix_Platy_TED_FBD_datation_3clock_40M.nex"
fichier_sortie="$DIR/../Analysis/Phylogeny/MrBayes/TED/Genera_phylogeny/Matrix_Platy_TED_FBD_datation_new_matrix_3clock_40M.nex "

regex=$(paste -sd'|' "$liste")
grep -E "$regex" "$fichier" > "$fichier_sortie"

echo "Extraction terminée : $fichier_sortie"
