################################################################################
# Name: helper_functions.R
# Authors: Lucas Buffan and Frédéric Lozes-Méléro
# Emails: lucas.l.buffan@gmail.com, frederic.lozes-melero@umontpellier.fr 
# Goal: Assistant function throughout our pipeline
################################################################################

library(seqinr)

## Standardised table writing --------------------------------------------------
write.tbl.std <- function(...){
  write.table(sep = "\t",
              na = "",
              row.names = FALSE,
              quote = FALSE,
              ...)
}

## Cut a string between two position (`from` and `to`) -------------------------
cut_string <- function(string, from, to){
  spl <- strsplit(string, split = "")[[1]]
  spl_before <- paste(spl[1:from], collapse = "")
  spl_after <- paste(spl[to:length(spl)], collapse = "")
  return(paste0(spl_before, spl_after))
}

## Convert seqinr alignement object to fasta -----------------------------------
fasta_from_aln <- function(Aln){
  aln_fasta <- ""
  for(i in 1:length(Aln$nam)){
    aln_fasta <- paste0(aln_fasta, ">", Aln$nam[i], "\n")
    aln_fasta <- paste0(aln_fasta, Aln$seq[[i]], "\n")
  }
  return(aln_fasta)
}

