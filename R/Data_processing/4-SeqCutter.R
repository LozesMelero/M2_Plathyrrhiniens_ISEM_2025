################################################################################
# Name: 4-SeqCutter.R
# Author: Frédéric Lozes-Méléro & Lucas Buffan
# Email: frederic.lozes-melero@umontpellier.fr & lucas.l.buffan@gmail.com 
# Goal: Cut alignement based on two indexed positions
#       /!\ USE IN COMMAND LINE /!\
################################################################################

library(seqinr)
library(DDD)

source("../useful/helper_functions.R")

## Specify arguments -----------------------------------------------------------
args <- commandArgs(trailingOnly = T)

path_to_aln <- args[1]       # Path to alignemnt
start <- as.numeric(args[2]) # Beginning of the portion of seq we want to remove
end <- as.numeric(args[3])   # Ending of the portion of seq we want to remove
out <- args[4]               # Path towards where we want to stored trimmed output

## Open alignment --------------------------------------------------------------
aln <- read.alignment(path_to_aln, format = "fasta")

## Cut sequence from position `start` to `end` ---------------------------------
aln_cut <- aln
aln_cut$seq <- lapply(X = aln_cut$seq,
                      FUN = cut_string,
                      from = start,
                      to = end)

## Save sequence after cutting -------------------------------------------------
aln_cut_fasta <- fasta_from_aln(aln_cut)
writeLines(aln_cut_fasta, con = out)
