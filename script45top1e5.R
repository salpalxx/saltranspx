#make qq-plot for each tissue
library(dplyr)
library(data.table)
library(ggplot2)
"%&%" = function(a,b) paste(a,b,sep="")
date <- Sys.Date()

########################
# To get the number of tested gene pairs, see how many lines are in the full output file and subtract 1
# from the command line:

# zcat GEUobs_v_GTExWB.meqtl.trans.diffchr.allres.txt.gz | wc -l
# 34389325

# Make a file of just the top 1e5 pairs for plotting
# zcat GEUobs_v_GTExWB.meqtl.trans.diffchr.allres.txt.gz | head -n 100001 > top1e5_GEUobs_v_GTExWB.meqtl.trans.diffchr.allres.txt


my.dir <- '/Users/sallyploch/mount/sally/GEUmeqtl/'
tis <- 'GEUobs_v_GTExWB'
resfile <- my.dir %&% 'top1e5_' %&% tis %&% '.meqtl.trans.diffchr.allres.txt.gz'
wcfile <- 34389325


res <- fread(resfile)

nn <- wcfile - 1
xx =  -log10((1:nn)/(nn+1)) #calculate expected based on the total tested

toplot <- dplyr::select(res,pvalue) %>% 
  mutate(obsP=sort(ifelse(pvalue<1e-30,30,-log10(pvalue)),decreasing=TRUE), expP=xx[1:1e5]) %>%
  dplyr::select(-pvalue)

png(filename = my.dir %&% "QQ_" %&% tis %&% ".png",res=100)
ggplot(toplot,aes(x=expP,y=obsP)) + geom_point(alpha=0.4) + coord_cartesian(xlim=c(-0.05,8.05),ylim=c(-0.05,30.05)) +
  geom_abline(slope=1,intercept = 0) + xlab(expression(Expected~~-log[10](italic(p)))) + 
  ylab(expression(Observed~~-log[10](italic(p)))) + ggtitle(tis)
dev.off()