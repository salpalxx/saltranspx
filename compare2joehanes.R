#read in diff chr results
library(dplyr)
library(data.table)
"%&%" = function(a,b) paste(a,b,sep="")

cohort <- "FHS"

tislist <- scan("/home/sally/tislist2.txt", "c")
jres <- fread("/home/sally/johanfiltered.csv", sep=",")
for(tis in tislist)
  {### load FHSobs v. GTExWBpred results P < 0.05
  fres <- fread('/home/sally/FHStranspx-cor_new/trans-diffChrs_FHSobs_v_GTEx' %&% tis %&% 'pred_pval0.05_2017-09-27.txt')
#print the number of trans-px gene pairs with P < 0.05
#print(dim(fres)[1])
#print the number of trans-px gene pairs with FDR < 0.05
fres_fdr0.05 <- dplyr::filter(fres,FDR<0.05)
transpairs <- dim(fres_fdr0.05)[1]

jres$SNP_Chr<-as.numeric(jres$SNP_Chr)
jres$Transcript_Chr<-as.numeric(jres$Transcript_Chr)

### compare Westra et al. eQTLs to trans-PrediXcan results
#join if chromosomes match
chrmatch <- inner_join(jres,fres,c("SNP_Chr"="pred_chr","Transcript_Chr"="obs_chr"))
#filter by position
posmatch <- dplyr::filter(chrmatch,SNP_Pos_hg19>pred_s1-1e6 & SNP_Pos_hg19 < pred_s2+1e6 & Transcript_Start_hg19 > obs_s1-1e6 & Transcript_Start_hg19 < obs_s2+1e6)
#filter by target genes matching
genematch <- dplyr::filter(posmatch,Transcript_GeneSymbol == obs_gene_id) %>% arrange(FDR)
#get number of Westra eQTLs that match trans-px gene pairs
#print(dim(genematch)[1])

#filter by unique trans-px gene pairs
unique_genematch <- genematch[!(duplicated(genematch$pred_ensembl) & duplicated(genematch$obs_ensembl)),]
#print(dim(unique_genematch)[1])

#get number of unique trans-px gene pairs with FDR<0.05
replicated <- dim(dplyr::filter(unique_genematch,FDR < 0.05))[1]

topmatches <- dplyr::filter(unique_genematch,FDR<0.05) %>% mutate(tis="FHSobs_v_GTEx" %&% tis %&% "pred")
#print(topmatches)
write.table(topmatches, "/home/sally/Joehanes_rep_trans-diffChrs_FHSobs_v_GTEx" %&% tis %&% "pred_FDR0.05_",quote=F,row.names=F)

results <- list(cohort, tis, transpairs, replicated)

#, trans-pairs, overlap)
fwrite(results, file="/home/sally/JoehanesRepTable", append=TRUE, sep=" ")

}

