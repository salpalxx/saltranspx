WB <- read.table("/Users/sallyploch/mount/sally/compare_cohorts/FDR0.05_FHSobs_v_GTExWB.meqtl.trans.diffchr_in_DGN_GEU.txt", header=TRUE)
repdgngeu <- select(WB, DGN_pval, GEU_pval)
repdgngeu %>% arrange(DGN_pval) #gives denominator for dgn 
repdgngeu %>% arrange(GEU_pval) #gives denominator for geu
filter(repdgngeu, DGN_pval < 0.05) #gives numerator for dgn
filter(repdgngeu, GEU_pval < 0.05) #gives numerator for geu






