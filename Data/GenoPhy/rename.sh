#!/bin/bash
# Script for renaming multiple taxa in FASTA files

mkdir -p genophy_cleaned_renamed  # create folder if it doesn't exist

for file in genophy_cleaned/*.fasta; do
    # Remplacer plusieurs noms de taxa
    sed -e 's/Alouatta_guariba_clamitans/Alouatta_guariba/g' \
        -e 's/Alouatta_seniculus_puruensis/Alouatta_puruensis/g' \
        -e 's/Ateles_fusciceps_rufiventris/Ateles_rufiventris/g' \
        -e 's/Cacajao_sp._a_FES-2022/Cacajao_amuna/g' \
        -e 's/Cacajao_calvus_novaesi/Cacajao_novaesi/g' \
        -e 's/Cacajao_calvus_rubicundus/Cacajao_rubicundus/g' \
        -e 's/Cacajao_calvus_ucayalii/Cacajao_ucayalii/g' \
        -e 's/Mico_humilis/Callibella_humilis/g' \
        -e 's/Chiropotes_israelita/Chiropotes_chiropotes/g' \
        -e 's/Leontocebus_fuscicollis_illigeri/Leontocebus_illigeri/g' \
        -e 's/Leontocebus_fuscicollis_lagonotus/Leontocebus_lagonotus/g' \
        -e 's/Leontocebus_fuscicollis_leucogenys/Leontocebus_leucogenys/g' \
        -e 's/Saguinus_melanoleucus_melanoleucus/Leontocebus_melanoleucus/g' \
        -e 's/Leontocebus_fuscicollis_nigrifrons/Leontocebus_nigrifrons/g' \
        -e 's/Leontocebus_fuscicollis_weddelli/Leontocebus_weddelli/g' \
        -e 's/Saguinus_geoffroyi/Oedipomidas_geoffroyi/g' \
        -e 's/Saguinus_leucopus/Oedipomidas_leucopus/g' \
        -e 's/Saguinus_oedipus/Oedipomidas_oedipus/g' \
        -e 's/Pithecia_pissinatti/Pithecia_pissinattii/g' \
        -e 's/Aotus_azarai/Aotus_azarae/g' \
        -e 's/Callicebus_caquetensis/Plecturocebus_caquetensis/g' \
        -e 's/Callicebus_ornatus/Plecturocebus_ornatus/g' \
        -e 's/Saimiri_sciureus_albigena/Saimiri_albigena/g' \
        -e 's/Saimiri_boliviensis_boliviensis/Saimiri_boliviensis/g' \
        -e 's/Sapajus_apella_macrocephalus/Sapajus_macrocephalus/g' \
        -e 's/Sapajus_nigritus_robustus/Sapajus_robustus/g' \
        -e 's/Saguinus_imperator/Tamarinus_imperator/g' \
        -e 's/Saguinus_inustus/Tamarinus_inustus/g' \
        -e 's/Saguinus_labiatus/Tamarinus_labiatus/g' \
        -e 's/Saguinus_labiatus_rufiventer/Tamarinus_labiatus/g' \
        -e 's/Saguinus_mystax/Tamarinus_mystax/g' \
        -e 's/Saguinus_imperator_subgrisescens/Tamarinus_subgrisescens/g' \
        -e 's/Callicebus_sp._n._GD-2019/Cheracebus_aquinoi/g' \
        "$file" > "genophy_cleaned_renamed/$(basename "$file" .fasta)_renamed.fasta"
done

echo "Task is completed"
