#Diffbind (Stark et al., 2011)

BiocManager::install("DiffBind")

library(DiffBind)

library(tidyverse)

setwd("/Users/ankitverma/Documents/Archivio2/tutorial/diffbind")


samples <- read.csv('short_sheet.csv')

names(samples)

#Create diffbind object

dbObj <- dba(sampleSheet=samples)

dbObj

#compute count information for each of the peaks/regions 

dbObj <- dba.count(dbObj, bUseSummarizeOverlaps=TRUE)

#Normalize

dbObj <- dba.normalize(dbObj)

#Establish contrast

dbObj <- dba.contrast(dbObj, categories=DBA_FACTOR, minMembers = 2)

#Performing the differential enrichment analysis

dbObj <- dba.analyze(dbObj, method=DBA_ALL_METHODS)

#Export DE sites

dbObj.DB <- dba.report(dbObj, method=DBA_DESEQ2, contrast = 1, th=1)

#Count Enriched and Depleted sites

sum(dbObj.DB$Fold>0)

sum(dbObj.DB$Fold<0)

#Create bed files for each keeping only significant peaks (p < 0.05)

out <- as.data.frame(dbObj.DB)

WT_enrich <- out %>% filter(FDR < 0.05 & Fold > 0) %>%  select(seqnames, start, end)

write.table(WT_enrich, file="WT_enriched.bed", sep="\t", quote=F, row.names=F, col.names=F)

KO_enrich <- out %>%  filter(FDR < 0.05 & Fold < 0) %>% select(seqnames, start, end)

write.table(KO_enrich, file="KO_enriched.bed", sep="\t", quote=F, row.names=F, col.names=F)

#Plot PCA

dba.plotPCA(dbObj,  attributes=DBA_FACTOR, label=DBA_ID, vColors = c("green","red"))

#Correlation heatmap

plot(dbObj)

#Assess

dba.show(dbObj, bContrasts=T)

dba.plotPCA(dbObj, contrast=1, method=DBA_DESEQ2, attributes=DBA_FACTOR, label=DBA_ID, vColors = c("red","green"))

#Plot venn

dba.plotVenn(dbObj,contrast=1,method=DBA_ALL_METHODS)

#Plot MA

dba.plotMA(dbObj, method=DBA_DESEQ2)

#Plot Volcano

dba.plotVolcano(dbObj)

#Heatmap

hmap <- colorRampPalette(c("red", "black", "green"))(n = 13)

readscores <- dba.plotHeatmap(dbObj, contrast=1, correlations=FALSE, scale="row", colScheme = hmap)
