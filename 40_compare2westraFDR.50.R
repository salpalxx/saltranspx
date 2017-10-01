#read in diff chr results
library(dplyr)
library(data.table)
library(ggplot2)
"%&%" = function(a,b) paste(a,b,sep="")

  cohort <- "FHS"

tislist <- scan("/home/sally/tislist.txt", "c")
wres <- fread("/home/wheelerlab1/Data/trans-px/Westra_results/2012-12-21-TransEQTLsFDR0.5_hg19.txt",sep="\t")
for(tis in tislist)
    {### load FHSobs v. GTExWBpred results P < 0.05
    fres <- fread('/home/sally/FHStranspx-cor_new/trans-diffChrs_FHSobs_v_GTEx' %&% tis %&% 'pred_pval0.05_2017-09-27.txt')
  #print the number of trans-px gene pairs with P < 0.05
  #print(dim(fres)[1])
  #print the number of trans-px gene pairs with FDR < 0.50
  fres_fdr0.50 <- dplyr::filter(fres,FDR<0.50)
  transpairs <- dim(fres_fdr0.50)[1]
  
  ### compare Westra et al. eQTLs to trans-PrediXcan results
  #join if chromosomes match
  chrmatch <- inner_join(wres,fres,c("Westra.SNPchr"="pred_chr","Westra.GENEchr"="obs_chr"))
  #filter by position
  posmatch <- dplyr::filter(chrmatch,Westra.SNPpos>pred_s1-1e6 & Westra.SNPpos < pred_s2+1e6 & Westra.GENEpos > obs_s1-1e6 & Westra.GENEpos < obs_s2+1e6)
  #filter by target genes matching
  genematch <- dplyr::filter(posmatch,Westra.GENEname == obs_gene_id) %>% arrange(FDR) %>% mutate(Westra.log10P=ifelse(Westra.P<1e-30,30,-log10(Westra.P)),transpx.log10P=ifelse(pvalue<1e-30,30,-log10(pvalue)))
  #get number of Westra eQTLs that match trans-px gene pairs
  #print(dim(genematch)[1])
  
  #filter by unique trans-px gene pairs
  unique_genematch <- genematch[!(duplicated(genematch$pred_ensembl) & duplicated(genematch$obs_ensembl)),]
  #print(dim(unique_genematch)[1])
  
  #get number of unique trans-px gene pairs with FDR<0.50
  replicated <- dim(dplyr::filter(unique_genematch,FDR < 0.50))[1]
  
    ggplot(genematch,aes(x=Westra.log10P,y=transpx.log10P)) + geom_point(alpha=0.4) + xlab("Westra et al. trans-eQTL -log10(P)") + ylab("FHS (GTEx-" %&% tis %&% ") trans-prediXcan -log10(P)") + geom_smooth(method='lm')
  summary(lm(data=genematch,transpx.log10P~Westra.log10P))
  
    topmatches <- dplyr::filter(unique_genematch,FDR<0.50) %>% mutate(tis="FHSobs_v_GTEx" %&% tis %&% "pred")
  #print(topmatches)
    write.table(topmatches, "/home/sally/Westra_rep_trans-diffChrs_FHSobs_v_GTEx" %&% tis %&% "pred_FDR0.50_",quote=F,row.names=F)
  
  results <- list(cohort, tis, transpairs, replicated)
  
  #, trans-pairs, overlap)
  fwrite(results, file="/home/sally/newFDR.50table.txt", append=TRUE, sep=" ")
  
}


