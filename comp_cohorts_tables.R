comp.50 <- comcohortstable
comp.50$RepInDGN..05 <- NULL
comp.50$RepinDGN..10 <- NULL
comp.50$RepinDGN..20 <- NULL
comp.50$RepInGEU..05 <- NULL
comp.50$RepinGEU..10 <- NULL
comp.50$RepinGEU..20 <- NULL
comp.50$DGNNum <- numdenoms$DGN.50Num
comp.50$DGNDenom <- numdenoms$DGN.50Denom
comp.50$GEUNum <- numdenoms$GEU.50Num
comp.50$GEUDenom <- numdenoms$GEU.50Denom
comp.50$DGNProp <- comp.50$DGNNum / comp.50$DGNDenom
comp.50$GEUProp <- comp.50$GEUNum / comp.50$GEUDenom
View(comp.50)
comp.50 <- comp.50[c(1,2,4,5,8,3,6,7,9)]
View(comp.50)
table.50 <- read.table("/Users/sallyploch/mount/sally/FDR.50.txt", header=TRUE)
table.50 <- filter(table.50, Cohort == "FHS")
table.50$Cohort <- NULL
table.50$DGNFrac <- comp.50$RepinDGN..50
table.50$DGNNum <- comp.50$DGNNum
table.50$DGNDenom <- comp.50$DGNNDenom
table.50$DGNProp <- comp.50$DGNNProp
table.50$GEUFrac <- comp.50$RepinGEU..50
table.50$GEUNum <- comp.50$GEUNum
table.50$GEUDenom <- comp.50$GEUDenom
table.50$GEUProp <- comp.50$GEUProp
table.50$GEUFrac <- NULL
table.50$GEUNum <- NULL
table.50$GEUDenom <- NULL
table.50$GEUProp <- NULL
table.50$DGNDenom <- comp.50$DGNDenom
table.50$DGNProp <- comp.50$DGNProp
table.50$GEUFrac <- comp.50$RepinGEU..50
table.50$GEUNum <- comp.50$GEUNum
table.50$GEUDenom <- comp.50$GEUDenom
table.50$GEUProp <- comp.50$GEUProp
write.table(table.50, file="/Users/sallyploch/mount/sally/compcohorts.50.txt", row.names = FALSE, col.names = TRUE, quote = FALSE)
