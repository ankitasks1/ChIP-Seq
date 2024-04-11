# Analysis pipeline for ChIP-Seq

<p><b> Do you have your own data? </b></p>

<p>Yes -> Proceed from Step 3  </p>
<p>No -> Follow from Step 1</p>


# Step 1: Explore Practice data, copy SRR code

https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE31477

eg.

GSM817343

USC_ChipSeq_HepG2_Input_UCDavis

GSM782122

USC_ChipSeq_HepG2_TCF7L2_UCDavis



# Step 2: Get the data from public repository (SRA)

#Install SRA tool kit

#https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit

#Select Ubuntu Linux 64 bit architecture/ mac OSX

#Install sra-toolkit

<code>/home/ankits/sratoolkit/bin/fastq-dump --split-files --gzip SRR9876543</code>


# Step 3:Quality Check your fastq files

#Get FastQC

#https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

#Click Download Now

#Click anyone as per your Operating system: FastQC v0.11.9 (Win/Linux zip file), FastQC v0.11.9 (Mac DMG image)

#Make sure the suitable java runtime environment (JRE) is installed : https://www.java.com/en/download/manual.jsp

<code>java -version</code>

#Install fastqc

#Make the fastqc function executable 

<code>chmod 755 fastqc</code>

<code>./fastqc SRR67548.fastq.gz</code>

#Expected output==> SRR67548.html


# Step 4:Trimming low quality and adapter sequences

#Trimmomatic

#Download 

<code>wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip</code>

<code>unzip Trimmomatic-0.39.zip</code>

#Single-End

<code>java -jar /pathTo/trimmomatic-0.35.jar SE -phred33 input.fq.gz output.fq.gz ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36</code>

#Paired-End

<code>java -jar /pathTo/trimmomatic-0.39.jar PE input_R1.fq.gz input_R2.fq.gz R1_paired.fastq.gz  R1_unpaired.fastq.gz  R2_paired.fastq.gz  R2_unpaired.fastq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36</code>


# Step 5: Mapping data: Alignment

#http://bowtie-bio.sourceforge.net/index.shtml

#Download the latest version: .zip

<code>unzip bowtie2-2.4.5-macos-arm64.zip</code>

<code>cd bowtie2-2.4.5-macos-arm64 </code>   

#Index the reference genome: 

<code>bowtie2-build GRCh38.fa GRCh38 </code>

#Alignment command:

#SE

<code>bowtie2 -q -x /home/ankits/folder/GRCh38 -U SRR2927818_trimmed.fastq -S SRR2927818.sam</code>

#PE

<code>bowtie2 -q -x home/ankits/folder/GRCh38 -1 SRR639251_paired.fastq -2 SRR639252_paired.fastq -U SRR639251_unpaired.fastq, SRR639252_unpaired.fastq -S SRR639252.sam</code>

# Step 6: Filter multimapping , duplicated and overlapping blacklisted region reads

#Samtools

http://www.htslib.org/download/ 

<code>cd samtools-1.xversion </code>

Install htslib first

<code>./configure --prefix=/where/to/install</code>

<code>make</code>

<code>make install</code>

Now install samtools similarly like htslib


#Convert SAM to BAM

<code>samtools view -S -b SRR639252.sam > SRR639252.bam</code>

#Extract Uniquely mapped reads

<code>samtools view -b -q 20 SRR639252.bam > SRR639252_uniq.bam</code>

#Sort BAM

<code>samtools sort -o SRR639252_uniq_sorted.bam SRR639252_uniq.bam</code>

<code>samtools index SRR639252_uniq_sorted.bam</code>

#Remove PCR duplicates

#https://github.com/broadinstitute/picard/releases/download/2.27.4/picard.jar

<code>java -jar ~/tools_av/picard.jar MarkDuplicates I=SRR639252_uniq_sorted.bam  O=SRR639252_uniq_sorted_dedups.bam M=SRR639252_picard_info.txt REMOVE_DUPLICATES=true AS=true CREATE_INDEX=true VALIDATION_STRINGENCY=STRICT</code>

#Removing “Blacklisted” regions

#Install Bedtools (required)

<code>wget https://github.com/arq5x/bedtools2/releases/download/v2.29.1/bedtools-2.29.1.tar.gz</code>

<code>tar -zxvf bedtools-2.29.1.tar.gz</code>

<code>cd bedtools2</code>

<code>make</code>

<code>bedtools intersect -abam SRR639252_uniq_sorted_dedups.bam -b hg38_blacklist.v2.bed -v | samtools sort --threads 4 -o SRR639252_uniq_sorted_dedups_freeblacklists.bam</code>

# Step 7: Peak calling

#Install MACS2 / SICER2

#How to get Conda

<code>wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh</code>

<code>./Miniconda3-latest-Linux-x86_64.sh</code>

<code>conda install -c bioconda macs2</code>

<code>conda install -c bioconda sicer2</code>

#Using Source

<code>wget https://github.com/macs3-project/MACS/archive/refs/tags/v2.2.7.1.tar.gz</code>

<code>tar zxvf MACS-2.2.6.tar.gz</code>

<code>cd MACS-2.2.6/</code>

<code>python setup.py install</code>

#Use SICER2 if you have histone marks 

https://files.pythonhosted.org/packages/6c/a8/25604c45b7eee1a56a0c459cd8b5edae9771a2669b60eb691aa16b055ee8/SICER2-1.0.3.tar.gz

<code>tar zxvf SICER2-1.0.3.tar.gz</code>

<code>cd SICER2-1.0.3</code>

<code>python setup.py install</code>

#Using pip

Require numpy

<code>pip install macs2</code>

<code>pip install SICER2</code>

#You need set path to run

Note: numpy in /Users/ankitverma/miniconda3/lib/python3.9/site-packages

Set path 

<code>export PATH=./miniconda3/:$PATH</code>

#Peak calling

#MACS2

<code>macs2 callpeak -t transcription_factor.bam -c Input.bam -g 2.7e9 -n transcription_factor --keep-dup all (add -f BAMPE for paired-end) </code>

#MACS2 --broad

<code>macs2 callpeak -t histone_marks.bam -c Input.bam -g 2.7e9 -n histone_marks --keep-dup all --broad (-f BAMPE for paired-end) </code>

#epic2

<code>epic2 -t histone_marks.bam -c Input.bam  --guess-bampe  --genome hg38 --output histone_marks --keep-duplicates</code>

#SICER2

<code>sicer --t histone_marks.bam -c Input.bam -s hg38 -w 200 -egf 0.85 -g 600 -o histone_marks</code>

# Step 8: Explore peaks

#Distribution of Peaks around features

#Deeptools

<code>computeMatrix reference-point -S H3K4me2.bw -R refTSS_v3.3_human_coordinate.hg38.bed -o test1 --a 3000 -b 3000 -bs 25 --missingDataAsZero</code>

<code>plotHeatmap -m test1 --colorList "white,blue" -out test1_computeMatrix1.png --sortUsing max </code>

#Motif analysis

#Web version https://meme-suite.org/meme/tools/meme

#Command-line version

#Install required packages, otherwise you will get error

Install XML perl module

<code>sudo apt-get install libxml-simple-perl</code>

Install HTML template perl module

<code>sudo apt-get install libhtml-template-perl</code>

Install JSON perl module

<code>sudo apt-get install libjson-perl</code>

Configure

<code>./configure --prefix=$HOME/meme --with-url=http://meme-suite.org --enable-build-libxml2 --enable-build-libxslt</code>

Cleanup old installations of meme

<code>make clean</code>

Install meme

<code>make</code>

<code>make test</code>

<code>make install</code>

#Get fasta (Obtain +/- 5bp from peak summit (generated by caller), you can use awk command to do that)

<code>bedtools getfasta -fi mm10.fa -bed peaks_summit_5bp.bed -fo out.fa</code>

eg. pf output fasta
Format of output: out.fa
>chr4:141410873-141410884
TGATTGATGCCTGCCTGTTA

#Run meme

<code>/home//meme/bin/meme-chip -meme-nmotifs 2 out.fa</code>

#Predict chromatin state

<code>wget  http://compbio.mit.edu/ChromHMM/ChromHMM.zip</code>

<code>java  -mx4000M -jar ~/tools_av/ChromHMM/ChromHMM.jar BinarizeBam  ~/tools_av/ChromHMM/CHROMSIZES/hg38.txt /bams chromhmmbamfiles.txt outputBinarizeddata</code>

<code>java  -mx4000M -jar ~/tools_av/ChromHMM/ChromHMM.jar LearnModel -p 2 outputBinarizeddata outLearnModel 10 hg38</code>

# Step 9: Annotate Peaks
#USE Linux based softwares

#####HOMER (http://homer.ucsd.edu/homer/introduction/install.html)

<code>annotatePeaks.pl peaks.bed hg38 > annotated_peaks.bed</code>
#or 
#USE R based softwares

#Here first you need to install R and RStudio

#R:https://cran.r-project.org/bin/windows/base/R-4.2.1-win.exe, https://cran.r-project.org/mirrors.html, https://cran.r-project.org/bin/macosx/base/R-4.2.1.pkg

#RStudio: https://www.rstudio.com/products/rstudio/download/

#Learn R scripting

#### Install Chipseeker
<pre>
```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("ChIPseeker")
BiocManager::install("ChIPseeker")
```
</pre>

#### Run Chipseeker	 
<pre>
```r
#Import library

library(ChIPseeker)

library(TxDb.Mmusculus.UCSC.mm9.knownGene)

library(clusterProfiler)

txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene

#Import Peak files 

SRR2927818_bed <- "/home/ankits/data_av/test/SRR2927818/peak_calling_bed/peak_annotation/chipseeker/peak_SRR2927817_818_peaks.bed"

SRR2927819_bed <- "/home/ankits/data_av/test/SRR2927818/peak_calling_bed/peak_annotation/chipseeker/peak_SRR2927817_819_peaks.bed"

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
```
</pre>

# Step 10: Differential binding analysis

<pre>
```r
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

</pre>

#### Diffbind 

<pre> 
ss <- read.csv("/mnt/beegfs6/home3/reid/av638/atacseq/surani_lab/wg/samplesheet.csv")
obj <- dba(sampleSheet=ss,  scoreCol=5, minOverlap=1)
# head(dba.peakset(obj, bRetrieve=TRUE, DataType=DBA_DATA_FRAME))
# dim(dba.peakset(obj, bRetrieve=TRUE, DataType=DBA_DATA_FRAME))
# colSums(dba.peakset(obj, minOverlap=2,bRetrieve=TRUE, DataType=DBA_DATA_FRAME)[,-c(1:3)])

# get counts (mapQCth=15 was kept to match with results from DESeq2/featurecounts, summits=FALSE for counting over full peak, )
obj <- dba.count(obj, minOverlap=1, score=DBA_SCORE_READS, summits=FALSE, mapQCth=15)# mapQC kept 15 for default of diffbind
obj <- dba.normalize(obj, method=DBA_DESEQ2, library = DBA_LIBSIZE_PEAKREADS) # peakreads # DBA_LIBSIZE_PEAKREADS for matching results the way DESEQ2 calculate size factor using library size withing the interval
obj <- dba.contrast(obj, minMembers = 2,categories=DBA_CONDITION) # create contrast using conditions present, atleast 2 replicates should be there for each condition
obj <- dba.analyze(obj) # run in default state

obj$norm$DESeq2$norm.calc
obj$norm$DESeq2$norm.method
obj$norm$DESeq2$lib.calc
obj$norm$DESeq2$lib.method
obj$norm$DESeq2$lib.sizes
obj$norm$DESeq2$norm.facs

overlapplot <- dba.overlap(obj, mode=DBA_OLAP_RATE)

dba.plotPCA(obj, attributes=DBA_CONDITION, label=DBA_REPLICATE, method= DBA_DESEQ2, score=DBA_SCORE_READS)

dba.plotPCA(obj, attributes=DBA_CONDITION, label=DBA_REPLICATE, method= DBA_DESEQ2, score=DBA_SCORE_NORMALIZED)

# normalized pca
db_pca_resn <- prcomp(t(counts(obj$DESeq2$DEdata, normalized =T)), center = TRUE, scale. = TRUE)
factoextra::fviz_pca_ind(db_pca_resn, geom = c("point", "text"), repel = TRUE, labelsize = 4)


df_diffbind_test_all <- data.frame(dba.report(obj, method=DBA_DESEQ2, contrast = 1, th=1))
df_diffbind_test_all["id"] <- paste0(df_diffbind_test_all$seqnames,
                                    "%",
                                    df_diffbind_test_all$start,
                                    "%",
                                    df_diffbind_test_all$end)
df_diffbind_test_all["antilogConc"] <- 2^df_diffbind_test_all$Conc
# df_diffbind_test_de <- df_diffbind_test_all %>% dplyr::filter(FDR < 0.05 & (Fold > log(2,2) | Fold < -log(2,2)))
df_diffbind_test_all0.05 <- df_diffbind_test_all %>% dplyr::filter(FDR < 0.05)
df_diffbind_test_de <- df_diffbind_test_all0.05 %>% dplyr::filter((Fold > 2 | Fold < -2))

df_diffbind_test_de["id"] <- paste0(df_diffbind_test_de$seqnames,
                                    "%",
                                    df_diffbind_test_de$start,
                                    "%",
                                    df_diffbind_test_de$end)


dim(df_diffbind_test_de)

test_contrasts <- dba.show(obj, bContrasts=TRUE)


</pre>
#### Deseq2 using Diffbind derived counts
<pre>

# get counts from diffbind
d_counts <- dba.peakset(obj, bRetrieve=TRUE, DataType=DBA_DATA_FRAME)
rownames(d_counts) <- paste0(d_counts$CHR, "%", d_counts$START, "%", d_counts$END)
d_counts <- d_counts[,-c(1:3)]
coldata <- data.frame(fread("coldata.txt", header = T))
rownames(coldata) <- coldata$SampleID
coldata$condition <- factor(coldata$Condition)
coldata$replicate <- factor(coldata$Replicate)
coldata <- coldata[,c(1,2,5,6)]
all(rownames(coldata) == colnames(counts)) #should print TRUE
# no filtering is required for low counts as the raw counts are already filtered by diffbind


d_dds <- DESeq2::DESeqDataSetFromMatrix(countData = as.matrix(d_counts), colData = coldata, design = ~ condition)

# no filtering is required for low counts as the raw counts are already filtered by diffbind
sizeFactors(d_dds) <- obj$DESeq2$facs
#-----
#test for local fittype
# d_dds <- estimateSizeFactors(d_dds) # alternative to DESeq2::DESeq(d_dds) for diffbind
d_dds <- estimateDispersions(d_dds, fitType = "local")
# DESeq2::nbinomWaldTest, with defaults
d_dds <- nbinomWaldTest(d_dds)

#-----
# vst pca
# d_vst <- DESeq2::vst(d_dds, fitType = "local")
# d_pcaData <- DESeq2::plotPCA(d_vst, intgroup=c("condition", "replicate"), returnData=TRUE)
# d_percentVar <- round(100 * attr(d_pcaData, "percentVar"))
# ggplot(d_pcaData, aes(PC1, PC2, color=condition, shape=replicate)) +
#   geom_point(size=3) +
#   xlab(paste0("PC1: ",d_percentVar[1],"% variance")) +
#   ylab(paste0("PC2: ",d_percentVar[2],"% variance")) +
#   coord_fixed() + theme_bw()
# 
# d_dds <- DESeq2::DESeq(d_dds)

# raw pca
d_pca_res <- prcomp(t(counts(d_dds)), center = TRUE, scale. = TRUE)
factoextra::fviz_pca_ind(d_pca_res, geom = c("point", "text"), repel = TRUE, labelsize =4)


# normalized pca
d_pca_resn <- prcomp(t(counts(d_dds, normalized =T)), center = TRUE, scale. = TRUE)
factoextra::fviz_pca_ind(d_pca_resn, geom = c("point", "text"), repel = TRUE, labelsize = 4)


res_d <- DESeq2::results(d_dds, contrast=c("condition", "DE", "hESC"))
res_d <- lfcShrink(d_dds, contrast=c("condition", "DE", "hESC"), res=res, type="ashr") # shrinkage as diffbind does for logfoldchange
res_d$threshold <- as.logical(res_d$padj < 0.05)
res_d_fdr <- data.frame(res_d[which(res_d$threshold == TRUE),])
# logfc filter
res_d_de <- res_d_fdr %>% dplyr::filter((log2FoldChange > 2) | (log2FoldChange < -2))
dim(res_d_de)
dim(df_diffbind_test_de)
length(intersect(rownames(res_d_de), df_diffbind_test_de$id))
res_d["id"] <- rownames(res_d)
res_d["logbaseMean"] <- log2(res_d$baseMean+1)
</pre>

#######----------Footnotes----------#### 


> obj <- dba.normalize(obj, method=DBA_DESEQ2)
> 
> obj$norm$DESeq2$norm.facs
> 
 [1] 1.1328567 1.1262136 1.3109969 1.1795663 1.3231150 1.0951867 1.1827066 1.1403031 0.7173917 0.6026628 0.7300126 0.6749297 0.8448875 0.8452929 1.2518313 0.7863995 1.3709245
[18] 1.1804729 0.7696433 0.7346065
> 
> obj <- dba.normalize(obj, method=DBA_DESEQ2,library = DBA_LIBSIZE_DEFAULT)
> 
> obj$norm$DESeq2$norm.facs
> 
 [1] 1.1328567 1.1262136 1.3109969 1.1795663 1.3231150 1.0951867 1.1827066 1.1403031 0.7173917 0.6026628 0.7300126 0.6749297 0.8448875 0.8452929 1.2518313 0.7863995 1.3709245
[18] 1.1804729 0.7696433 0.7346065
> 
> obj <- dba.normalize(obj, method=DBA_DESEQ2, normalize = DBA_NORM_DEFAULT)
> obj$norm$DESeq2$norm.facs
> 
 [1] 1.1328567 1.1262136 1.3109969 1.1795663 1.3231150 1.0951867 1.1827066 1.1403031 0.7173917 0.6026628 0.7300126 0.6749297 0.8448875 0.8452929 1.2518313 0.7863995 1.3709245
[18] 1.1804729 0.7696433 0.7346065

> obj <- dba.normalize(obj, library = DBA_LIBSIZE_FULL)
> 
> obj$norm$DESeq2$norm.facs
> 
 [1] 1.1328567 1.1262136 1.3109969 1.1795663 1.3231150 1.0951867 1.1827066 1.1403031 0.7173917 0.6026628 0.7300126 0.6749297 0.8448875 0.8452929 1.2518313 0.7863995 1.3709245
[18] 1.1804729 0.7696433 0.7346065

> obj <- dba.normalize(obj, method=DBA_DESEQ2, normalize = DBA_NORM_NATIVE,library = DBA_LIBSIZE_DEFAULT)
> 
> obj$norm$DESeq2$norm.facs
> 
        DE_r1         DE_r2       hESC_r1       hESC_r2  hPGCLC_d2_r1  hPGCLC_d2_r2  hPGCLC_d4_r1  hPGCLC_d4_r2 hPGC_wk7_F_r1 hPGC_wk7_F_r2 hPGC_wk7_F_r3 hPGC_wk7_M_r1 
    1.3996433     1.8616030     1.4199660     2.3250907     1.1248697     1.7218136     1.2156686     1.4979737     0.5409706     0.5401057     0.5342037     0.6621932 
hPGC_wk7_M_r2         ME_r1         ME_r2         ME_r3      PreME_r1      PreME_r2 Soma_wk7_M_r1 Soma_wk7_M_r2 
    0.7019794     0.6122090     1.5503055     1.1209880     1.3732644     1.6133297     0.4713268     0.6442923

> obj <- dba.normalize(obj, method=DBA_DESEQ2, normalize = DBA_NORM_NATIVE)# RLE for DESEq2
> 
> obj$norm$DESeq2$norm.facs
> 
        DE_r1         DE_r2       hESC_r1       hESC_r2  hPGCLC_d2_r1  hPGCLC_d2_r2  hPGCLC_d4_r1  hPGCLC_d4_r2 hPGC_wk7_F_r1 hPGC_wk7_F_r2 hPGC_wk7_F_r3 hPGC_wk7_M_r1 
    1.3996433     1.8616030     1.4199660     2.3250907     1.1248697     1.7218136     1.2156686     1.4979737     0.5409706     0.5401057     0.5342037     0.6621932 
hPGC_wk7_M_r2         ME_r1         ME_r2         ME_r3      PreME_r1      PreME_r2 Soma_wk7_M_r1 Soma_wk7_M_r2 
    0.7019794     0.6122090     1.5503055     1.1209880     1.3732644     1.6133297     0.4713268     0.6442923 

> obj <- dba.normalize(obj, method=DBA_DESEQ2, normalize = DBA_NORM_NATIVE,library = DBA_LIBSIZE_PEAKREADS)
> 
> obj$norm$DESeq2$norm.facs
> 
        DE_r1         DE_r2       hESC_r1       hESC_r2  hPGCLC_d2_r1  hPGCLC_d2_r2  hPGCLC_d4_r1  hPGCLC_d4_r2 hPGC_wk7_F_r1 hPGC_wk7_F_r2 hPGC_wk7_F_r3 hPGC_wk7_M_r1 
    1.3996433     1.8616030     1.4199660     2.3250907     1.1248697     1.7218136     1.2156686     1.4979737     0.5409706     0.5401057     0.5342037     0.6621932 
hPGC_wk7_M_r2         ME_r1         ME_r2         ME_r3      PreME_r1      PreME_r2 Soma_wk7_M_r1 Soma_wk7_M_r2 
    0.7019794     0.6122090     1.5503055     1.1209880     1.3732644     1.6133297     0.4713268     0.6442923 
> obj <- dba.normalize(obj, method=DBA_DESEQ2, normalize = DBA_NORM_DEFAULT,library = DBA_LIBSIZE_PEAKREADS)
> 
> obj$norm$DESeq2$norm.facs
> 
        DE_r1         DE_r2       hESC_r1       hESC_r2  hPGCLC_d2_r1  hPGCLC_d2_r2  hPGCLC_d4_r1  hPGCLC_d4_r2 hPGC_wk7_F_r1 hPGC_wk7_F_r2 hPGC_wk7_F_r3 hPGC_wk7_M_r1 
    1.3893406     1.9286069     1.1327701     1.9583004     0.9139431     1.6062935     1.0193394     1.2666328     0.4496935     0.4864651     0.4525491     0.6153672 
hPGC_wk7_M_r2         ME_r1         ME_r2         ME_r3      PreME_r1      PreME_r2 Soma_wk7_M_r1 Soma_wk7_M_r2 
    0.6240718     0.4666545     1.2619604     0.9095477     1.1036117     1.3171674     0.4056359     0.6920489
> 
> obj <- dba.normalize(obj, method=DBA_DESEQ2,library = DBA_LIBSIZE_PEAKREADS)
> 
> obj$norm$DESeq2$norm.facs
> 
        DE_r1         DE_r2       hESC_r1       hESC_r2  hPGCLC_d2_r1  hPGCLC_d2_r2  hPGCLC_d4_r1  hPGCLC_d4_r2 hPGC_wk7_F_r1 hPGC_wk7_F_r2 hPGC_wk7_F_r3 hPGC_wk7_M_r1 
    1.3893406     1.9286069     1.1327701     1.9583004     0.9139431     1.6062935     1.0193394     1.2666328     0.4496935     0.4864651     0.4525491     0.6153672 
hPGC_wk7_M_r2         ME_r1         ME_r2         ME_r3      PreME_r1      PreME_r2 Soma_wk7_M_r1 Soma_wk7_M_r2 
    0.6240718     0.4666545     1.2619604     0.9095477     1.1036117     1.3171674     0.4056359     0.6920489
> 
> obj <- dba.normalize(obj, method=DBA_DESEQ2, normalize = DBA_NORM_NATIVE,library = DBA_LIBSIZE_FULL)
> 
> obj$norm$DESeq2$norm.facs
> 
        DE_r1         DE_r2       hESC_r1       hESC_r2  hPGCLC_d2_r1  hPGCLC_d2_r2  hPGCLC_d4_r1  hPGCLC_d4_r2 hPGC_wk7_F_r1 hPGC_wk7_F_r2 hPGC_wk7_F_r3 hPGC_wk7_M_r1 
    1.3996433     1.8616030     1.4199660     2.3250907     1.1248697     1.7218136     1.2156686     1.4979737     0.5409706     0.5401057     0.5342037     0.6621932 
hPGC_wk7_M_r2         ME_r1         ME_r2         ME_r3      PreME_r1      PreME_r2 Soma_wk7_M_r1 Soma_wk7_M_r2 
    0.7019794     0.6122090     1.5503055     1.1209880     1.3732644     1.6133297     0.4713268     0.6442923
    
> obj <- dba.normalize(obj,normalize = DBA_NORM_NATIVE)
> 
> obj$norm$DESeq2$norm.facs
        DE_r1         DE_r2       hESC_r1       hESC_r2  hPGCLC_d2_r1  hPGCLC_d2_r2  hPGCLC_d4_r1  hPGCLC_d4_r2 hPGC_wk7_F_r1 hPGC_wk7_F_r2 hPGC_wk7_F_r3 hPGC_wk7_M_r1 
    1.3996433     1.8616030     1.4199660     2.3250907     1.1248697     1.7218136     1.2156686     1.4979737     0.5409706     0.5401057     0.5342037     0.6621932 
hPGC_wk7_M_r2         ME_r1         ME_r2         ME_r3      PreME_r1      PreME_r2 Soma_wk7_M_r1 Soma_wk7_M_r2 
    0.7019794     0.6122090     1.5503055     1.1209880     1.3732644     1.6133297     0.4713268     0.6442923
> 
> obj <- dba.normalize(obj)# default all
> 
> obj$norm$DESeq2$norm.facs
> 
 [1] 1.1328567 1.1262136 1.3109969 1.1795663 1.3231150 1.0951867 1.1827066 1.1403031 0.7173917 0.6026628 0.7300126 0.6749297 0.8448875 0.8452929 1.2518313 0.7863995 1.3709245
[18] 1.1804729 0.7696433 0.7346065

> obj <- dba.normalize(obj, library = DBA_LIBSIZE_PEAKREADS)
> 
> obj$norm$DESeq2$norm.facs
> 
        DE_r1         DE_r2       hESC_r1       hESC_r2  hPGCLC_d2_r1  hPGCLC_d2_r2  hPGCLC_d4_r1  hPGCLC_d4_r2 hPGC_wk7_F_r1 hPGC_wk7_F_r2 hPGC_wk7_F_r3 hPGC_wk7_M_r1 
    1.3893406     1.9286069     1.1327701     1.9583004     0.9139431     1.6062935     1.0193394     1.2666328     0.4496935     0.4864651     0.4525491     0.6153672 
hPGC_wk7_M_r2         ME_r1         ME_r2         ME_r3      PreME_r1      PreME_r2 Soma_wk7_M_r1 Soma_wk7_M_r2 
    0.6240718     0.4666545     1.2619604     0.9095477     1.1036117     1.3171674     0.4056359     0.6920489

> obj <- dba.normalize(obj, method=DBA_DESEQ2, library = DBA_LIBSIZE_FULL)
> 
> 
> obj$norm$DESeq2$norm.facs
> 
 [1] 1.1328567 1.1262136 1.3109969 1.1795663 1.3231150 1.0951867 1.1827066 1.1403031 0.7173917 0.6026628 0.7300126 0.6749297 0.8448875 0.8452929 1.2518313 0.7863995 1.3709245
[18] 1.1804729 0.7696433 0.7346065

> obj <- dba.normalize(obj, method=DBA_DESEQ2, library = DBA_LIBSIZE_PEAKREADS)
> 
> obj$norm$DESeq2$norm.facs
        DE_r1         DE_r2       hESC_r1       hESC_r2  hPGCLC_d2_r1  hPGCLC_d2_r2  hPGCLC_d4_r1  hPGCLC_d4_r2 hPGC_wk7_F_r1 hPGC_wk7_F_r2 hPGC_wk7_F_r3 hPGC_wk7_M_r1 
    1.3893406     1.9286069     1.1327701     1.9583004     0.9139431     1.6062935     1.0193394     1.2666328     0.4496935     0.4864651     0.4525491     0.6153672 
hPGC_wk7_M_r2         ME_r1         ME_r2         ME_r3      PreME_r1      PreME_r2 Soma_wk7_M_r1 Soma_wk7_M_r2 
    0.6240718     0.4666545     1.2619604     0.9095477     1.1036117     1.3171674     0.4056359     0.6920489


> obj <- dba.normalize(obj, method=DBA_DESEQ2, library = DBA_LIBSIZE_DEFAULT)
> 
> obj$norm$DESeq2$norm.facs
> 
[1] 1.1328567 1.1262136 1.3109969 1.1795663 1.3231150 1.0951867 1.1827066 1.1403031 0.7173917 0.6026628 0.7300126 0.6749297 0.8448875 0.8452929 1.2518313 0.7863995 1.3709245
[18] 1.1804729 0.7696433 0.7346065



