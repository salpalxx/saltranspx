library(dplyr)
library(data.table)
library(tibble)
"%&%" = function(a,b) paste(a,b,sep="")

#Volumes: osx, group: linux
predfile <- "/Users/sallyploch/mount/sally/FHSpredexp/WB_test_predicted_expression.txt"
obsfile <- "/Users/sallyploch/mount/wheelerlab1/Data/trans-px/obsFHSexp/FHS_obs_expression.txt.gz"
annotfile <- "/Users/sallyploch/mount/wheelerlab1/Data/trans-px/gencode.v18.genes.patched_contigs.summary.protein"

#gzcat: osx , gzcat: linux
pred <- fread('gzcat ' %&% predfile)
obs <- fread('gzcat ' %&% obsfile)
annot <- fread(annotfile)

#convert to matrix & assign rownames (IID) and colnames (ENSEMBL w/o version)
predmat <- as.matrix(pred[,-1:-2])
rownames(predmat) <- pred$IID
colnames(predmat) <- substr(colnames(predmat),1,15) #ENSG w/o version
obsmat <- as.matrix(obs[,-1])
rownames(obsmat) <- obs$IID

#get intersection vectors
samplelist <- as.character(intersect(pred$IID, obs$IID))
genelist <- as.character((intersect(colnames(predmat),colnames(obsmat))))

#get intersection matrices
predmat <- predmat[samplelist,genelist]
obsmat <- obsmat[samplelist,genelist]

#check matrices match
stopifnot(colnames(predmat) == colnames(obsmat))
stopifnot(rownames(predmat) == rownames(obsmat))

#make location files
annotdf <- as.data.frame(annot)
annotdf <- mutate(annotdf, gene = substr(V5,1,15)) %>% 
  dplyr::filter(gene %in% genelist) %>% arrange(gene) %>%
  dplyr::rename(chr=V1,s1=V3,s2=V4) %>% 
  dplyr::mutate(chrnum = as.numeric(substr(chr,4,5))) %>% arrange(chrnum,s1) #sort by chr, pos
obsloc <- dplyr::select(annotdf,gene,chr,s1,s2)
predloc <- dplyr::select(obsloc,gene,chr,s1)

#sort matrices by chr, pos
newordergenes <- obsloc$gene
predmat <- predmat[,newordergenes]
obsmat <- obsmat[,newordergenes]

#remove genes w/sum = 0 for pred exp
keepgenes = colSums(predmat) != 0
predmat <- predmat[,keepgenes]
obsmat <- obsmat[,keepgenes]
predloc <- predloc[keepgenes,]
obsloc <- obsloc[keepgenes,]

#check matrices match
stopifnot(colnames(predmat) == colnames(obsmat))
stopifnot(rownames(predmat) == rownames(obsmat))
stopifnot(colnames(predmat) == predloc$gene)
stopifnot(colnames(obsmat) == obsloc$gene)

#transpose matrices to MatrixEQTL format
tpred <- as.data.frame(t(predmat))
tpred <- rownames_to_column(tpred, var="gene")
tobs <- as.data.frame(t(obsmat))
tobs <- rownames_to_column(tobs, var="gene")

#make gene annot file
annot2write <- mutate(annotdf,ensembl_id = gene, gene_id = V6, strand = V2) %>% 
  dplyr::select(chr, strand, s1, s2, ensembl_id, gene_id, chrnum )

write.table(predloc, file="/Users/sallyploch/mount/sally/FHSmeqtlformat/predLOC_gtexWB.txt", quote=F,row.names = F,sep='\t')
write.table(obsloc, file="/Users/sallyploch/mount/sally/FHSmeqtlformat/obsLOC_gtexWB.txt", quote=F,row.names = F,sep='\t')
write.table(tpred, file="/Users/sallyploch/mount/sally/FHSmeqtlformat/predEXP_FHS_from_gtexWB.txt", quote=F,row.names = F,sep='\t')
write.table(tobs, file="/Users/sallyploch/mount/sally/FHSmeqtlformat/obsEXP_FHS_gtexWB.txt", quote=F,row.names = F,sep='\t')
write.table(annot2write, file="/Users/sallyploch/mount/sally/FHSmeqtlformat/FHS_gene_annot_gtexWB.txt", quote=F, row.names=F, sep='\t')