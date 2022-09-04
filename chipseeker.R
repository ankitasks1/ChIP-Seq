#-------Peak annotation----------
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ChIPseeker")
BiocManager::install("ChIPseeker")

#Run Chipseeker	 

#Import library

library(ChIPseeker)

library(TxDb.Mmusculus.UCSC.mm9.knownGene)

library(clusterProfiler)

txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene

#Import Peak files 

SRR2927818_bed <- "/home/ankits/data_av/test/SRR2927818/peak_calling_bed/peak_annotation/chipseeker/peak_SRR2927817_818_peaks.bed"

SRR2927819_bed <- "/home/ankits/data_av/test/SRR2927818/peak_calling_bed/peak_annotation/chipseeker/peak_SRR2927817_819_peaks.bed.gz"

print(SRR2927818_bed)

SRR2927818_peak <- readPeakFile(SRR2927818_bed)

SRR2927818_peak

SRR2927819_peak <- readPeakFile(SRR2927819_bed)

SRR2927819_peak

#ChIP peaks coverage plot

covplot(SRR2927819_peak, weightCol="V5")

covplot(SRR2927819_peak, weightCol="V5", chrs=c("chr17", "chr18"), xlim=c(4.5e7, 5e7))

#Profile of ChIP peaks binding to TSS regions

promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)

SRR2927819_tagMatrix <- getTagMatrix(SRR2927819_peak, windows=promoter)

#Heatmap of ChIP binding to TSS regions

tagHeatmap(SRR2927819_tagMatrix, xlim=c(-3000, 3000), color="red")

peakHeatmap(SRR2927819_bed, TxDb=txdb, upstream=3000, downstream=3000, color="red")

#Average Profile of ChIP peaks binding to TSS region

plotAvgProf(SRR2927819_tagMatrix, xlim=c(-3000, 3000),
            xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")

plotAvgProf2(SRR2927819_bed, TxDb=txdb, upstream=3000, downstream=3000,
             xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")

plotAvgProf(SRR2927819_tagMatrix, xlim=c(-3000, 3000), conf = 0.95, resample = 1000)


#Peak Annotation

SRR2927819_peakAnno <- annotatePeak(SRR2927819_bed,, tssRegion=c(-3000, 3000),
                                    TxDb=txdb, annoDb="org.Mm.eg.db")

#Visualize Genomic Annotation

plotAnnoPie(SRR2927819_peakAnno)

plotAnnoBar(SRR2927819_peakAnno)

vennpie(SRR2927819_peakAnno)

upsetplot(SRR2927819_peakAnno)

upsetplot(SRR2927819_peakAnno, vennpie=TRUE)

#Visualize distribution of TF-binding loci relative to TSS

plotDistToTSS(SRR2927819_peakAnno,
              title="Distribution of transcription factor-binding loci\nrelative to TSS")

