


library(tidyverse)


table <- read.delim("./Data/GenoPhy/accession.txt", header = F)

genes <- c("ABCA1", "ADORA3", "ERC2", "FES", "FOXP1", "MAPKAP1", 
           "RAG1", "RPGRIP1", "SIM1", "ZFX", "cOX1", "cOX2", "cOX3", "cYTB", "nD2")

f_table <- table[grepl(paste(genes, collapse = "|"), table$V2), ]

f_table$V2 <- sub(".*/", "", f_table$V2)  
f_table$V2 <- sub("dir$", "", f_table$V2)  
f_table$V1 <- sub(".*/", "", f_table$V1)

fwid_table <- f_table %>% pivot_wider(names_from = V2, values_from = V3)
fwid_table <- fwid_table[,1:16]

write.table(fwid_table, "../../../fichier_table.txt", sep = "\t", row.names = FALSE, quote = FALSE)

#### Cata
table <- read.delim("./Data/GenoPhy/accession_cata.txt", header = F)

genes <- c("ABCA1", "ADORA3", "ERC2", "FES", "FOXP1", "MAPKAP1", 
           "RAG1", "RPGRIP1", "SIM1", "ZFX", "cOX1", "cOX2", "cOX3", "cYTB", "nD2")

f_table <- table[grepl(paste(genes, collapse = "|"), table$V2), ]

f_table$V2 <- sub(".*/", "", f_table$V2)  
f_table$V2 <- sub("dir$", "", f_table$V2)  
f_table$V1 <- sub(".*/", "", f_table$V1)

fwid_table <- f_table %>% pivot_wider(names_from = V2, values_from = V3)
fwid_table <- fwid_table[,1:16]

write.table(fwid_table, "../../../fichier_table.txt", sep = "\t", row.names = FALSE, quote = FALSE)














