#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 19 12:50:21 2025

This script is meant to facilitate sequence re-indexing in the .nex file after 
reducing the morphological matrix.

@author: lucas.buffan
"""

## Useful functions
def paste(L, sep = " "):
    final_str = str()
    for el in L[:-1]: # avoids ending with a separator
        final_str += el + sep
    return(final_str+L[-1])

def decrease_doublet(doublet, by): # decrease a doublet of index like `123-140` by a given step
    # Lower bound    
    low = int(doublet.split("-")[0])
    if low > 1: #don't remove anything on the first position:
        low -= diff
    # Upper bound
    delim = ""
    up = doublet.split("-")[1]
    if up[-2:] == ";\n":
        delim = ";\n"
        up = up[:-2]
    elif up[-1:] == ";":
        delim = ";"
        up = up[:-1]
    up = int(up) - diff
    # Assemble
    return(str(low) + "-" + str(up) + delim)


## Ctype ordering

frame = 85 # frame shift length

# Former "c-typing" of morphological characters
A = "86 90 - 92 94 96 - 97 102 - 106 109 - 111 113 - 114 118 - 119 122 - 123 128 - 133 136 139 143 145 - 149 151 - 152 154 164 - 167 169 171 - 174 176 - 177 187 - 190 192 - 198 204 206 - 210 213 215 - 218 221 - 222 228 - 236 241 - 244 247 - 253 256 - 260 263 - 267 270 274 276 - 282 289 293 - 294 298 - 299 302 - 308 310 - 313 315 317 319 321 - 322 324 326 - 328 330 - 335 337 - 344 351 - 354 356 - 357 370"
B = A.split(" ")
C = []
for el in B:
    if el != '-':
        el = str(int(el) - frame)
    C.append(el)
    
fin_string = paste(C)


## Partition indices

diff = 168 # = 456 - 288
new_part = ""
with open("/home/lucas.buffan/Documents/GitHub/PhD/Chapter_3/Tambouille_Lucas/old_partitions.nex") as file:
    for line in file:
        if "charset" in line:
            spl_line = line.split(" = ")
            positions = spl_line[1]
            # If line only contains one position doublet
            if len(positions.split(" ")) == 1:
                new_part += spl_line[0] + " = " + decrease_doublet(positions, by = diff)
            # Otherwise, consider each doublet one by one
            else:
                pos_list = positions.split(" ")
                new_pos_list = [decrease_doublet(x, by = diff) for x in pos_list]
                new_part += spl_line[0] + " = " + paste(new_pos_list)
        else:
            new_part += line

# Save new partitions
New_partitions = open("/home/lucas.buffan/Documents/GitHub/PhD/Chapter_3/Tambouille_Lucas/new_partitions.nex", "w")
n = New_partitions.write(new_part)
New_partitions.close()