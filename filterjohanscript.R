library(data.table)
library(dplyr)
johan <- fread("/home/sally/eqtl-gene-1000g-peer-validate_c1-23_c28.csv")
johan <- filter(johan, SNP_Chr != Transcript_Chr)
johan <- filter(johan, log10FDR <= log10(0.05))
cols <- c(2,3,4,12,14,16,22,23)
johannew <- johan[,cols]
write.csv(johannew, file = "/home/sally/johanfiltered.csv", append = FALSE, quote = FALSE, row.names= FALSE, col.names = TRUE)