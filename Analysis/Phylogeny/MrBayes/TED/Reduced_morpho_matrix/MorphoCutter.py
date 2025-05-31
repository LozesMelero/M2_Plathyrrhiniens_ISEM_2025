#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 18 10:17:24 2025

Subdivide morphological character matrix (Laurent's Nightmare) based on indices

@author: lucas.buffan
"""

import os

os.chdir("/home/lucas.buffan/Documents/GitHub/PhD/Chapter_3/Platy_TED/Analysis/Phylogeny/MrBayes/TED/")

## GPT-flavoured fucntion to return the index of the character each element from a morphological sequence belongs to (including parentheses/brackets)
        # e.g. "01{01}0" => [0, 1, 2,2,2,2, 3]
def parse_morpho_seq(morpho_seq):
    char_indices = [] # will store the index of the character each state encoded in the sequence belongs to
    current_char = 0
    i = 0

    while i < len(morpho_seq):
        if morpho_seq[i] in '({':
            # set the closing statement according to the string element handeled
            closing = ')' if morpho_seq[i] == '(' else '}'
            start = i + 1
            end = morpho_seq.find(closing, start) # returns the index of the first occurrence of the closing statement starting from the opening one

            # Assign current_char to each state inside the parentheses/brackets + those of the parentheses/brackets
            num_states = end - start + 2
            char_indices.extend([current_char] * num_states)

            # Skip past the closing bracket
            i = end + 1
            current_char += 1
        else:
            # Single state character
            char_indices.append(current_char)
            i += 1
            current_char += 1

    return char_indices


# First and last positions of the characters we want to include (restrict to the 288 dental characters)
first = 85 # First dental character (nb. indexing starts at 0 !)
last = 372 # Last dental character

# Slight changes in the header
n_morpho_characters = last - first
tot_length = 12674-(456-n_morpho_characters)+1
l4 = "	DIMENSIONS NTAX = 250 NCHAR = " + str(tot_length) + ";\n"
l5 = "	FORMAT DATATYPE = mixed(standard:1-" + str(n_morpho_characters+1) + ", DNA:" + str(n_morpho_characters+2) + "-" + str(tot_length) + ") interleave=yes  GAP = - MISSING = ?;\n"

reduced_mat = str()
with open("Platy_TED_FBD_datation.nex") as nex:
    line_idx = 0
    for line in nex:
        # Make sure we restrict to morpho
        if "[MOLECULAR]" in line:
            break
        # Leave header as it is
        if line_idx < 10 or line == "\n":
            if "DIMENSIONS" in line:
                reduced_mat += l4
            elif "DATATYPE" in line:
                reduced_mat += l5
            else:
                reduced_mat += line
        # Work on morphological sequences
        else:
            line_spl = line.split(" ")
            tax_name = line_spl.pop(0) # first item = taxon name
            morpho_seq = [i for i in line_spl if len(i) > 0]
            # Retrieve the entire morphological sequence without separator tho
            full_morpho_seq = ''.join(morpho_seq)[:-1]
            
            # Assign each element a character number 
            final_indices = parse_morpho_seq(full_morpho_seq)
            
            # Get the position of the characters within the range we want to retain
            retained_indices = [i for i, x in enumerate(final_indices) if first <= x <= last]
            # Filter out non-dental characters. As it's a consecutive range, we just need the first and last elements of retained_indices
            dental_characters = full_morpho_seq[retained_indices[0]:(retained_indices[-1]+1)]
            # Add species name and store
            reduced_mat += tax_name + "\t" + dental_characters + "\n"
        line_idx += 1

# Add terminal semicolon
reduced_mat += ";"

# Save it
Dental_matrix = open("../../../../../Tambouille_Lucas/Dental_matrix_only.nex", "w")
n = Dental_matrix.write(reduced_mat)
Dental_matrix.close()