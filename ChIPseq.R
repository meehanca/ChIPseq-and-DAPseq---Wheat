library(ChIPpeakAnno)
library(ChIPseeker)
library(biomaRt)
library(rtracklayer)
library(GenomicFeatures)
library(GenomicDistributions)
library(GenomicDistributionsData)
library(GenomicRanges)
library(circlize)
library(data.table)
library(EnrichedHeatmap)
library(magick)
library(ComplexHeatmap)
library(dplyr)


listMarts(host="plants.ensembl.org")
mart <- useMart("plants_mart", host = "plants.ensembl.org")
listDatasets(mart)
ensembl <- useDataset(dataset = "taestivum_eg_gene", mart = mart)
txdb <- makeTxDbFromBiomart("plants_mart","taestivum_eg_gene",host="plants.ensembl.org")

####################################
### ChIP-seq ###
####################################

####################################
### MACS2 - RKD ###
####################################

setwd("~/Documents/ChIPseq")

file_narrowPeak = "MACS2/ChIPseq_RKD2q0.05_peaks.narrowPeak"

extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric",
                          qValue = "numeric", peak = "integer")

RKD_ChIP_MACS2_narrowPeak <- import(file_narrowPeak, format = "BED",
                           extraCols = extraCols_narrowPeak)

RKD_ChIP_Wheat_MACS2_annotation <- annotatePeakInBatch(RKD_ChIP_MACS2_narrowPeak,ensembl,"TSS", output = "nearestBiDirectionalPromoters",bindingRegion = c(-10000, 10000))

####################################
### MACS2 - RKD2 ###
####################################

setwd("~/Documents/ChIPseq")

file_narrowPeak = "MACS2/ChIPseq_controlq0.05_peaks.narrowPeak"
extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric",
                          qValue = "numeric", peak = "integer")

RKD2_ChIP_MACS2_narrowPeak <- import(file_narrowPeak, format = "BED",
                                    extraCols = extraCols_narrowPeak)

RKD2_ChIP_Wheat_MACS2_annotation <- annotatePeakInBatch(RKD2_ChIP_MACS2_narrowPeak,ensembl,"TSS", output = "nearestBiDirectionalPromoters",bindingRegion = c(-10000, 10000))

####################################
### MACS2 - Overlap ###
####################################
RKD2_ChIP1 <- RKD_ChIP_MACS2_narrowPeak
RKD2_ChIP2 <- RKD2_ChIP_MACS2_narrowPeak

Overlap <- findOverlapsOfPeaks(RKD1_ChIP,RKD2_ChIP)
makeVennDiagram(Overlap,Txdb=txdb,by="feature",)

ChIPseq_MACS2_narrowPeak <- makeGRangesFromDataFrame(Overlap$overlappingPeaks)

Annotation<-getAnnotation(ensembl)

RKD_ChIP_Wheat_MACS2_annotation <- annotatePeakInBatch(ChIPseq_MACS2_narrowPeak,ensembl,"TSS",Annotation, output = "nearestBiDirectionalPromoters",bindingRegion = c(-10000, 10000))

write.table(RKD_ChIP_Wheat_MACS2_annotation,"RKD_ChIP_MACS2_Wheat_peaks.txt",quote=F,sep="\t")
write.table(RKD_ChIP_Wheat_MACS2_annotation$feature,"RKD_ChIP_seq_MACS2_gene_IDs.txt",quote=F,row.names = F,col.names = F)

####################################
### DAPseq ###
####################################

####################################
### MACS2 ###
####################################

setwd("~/Documents/DAPseq/Wheat")

file_narrowPeak = "Run with correct genome size/RKD2_HALOq_peaks.narrowPeak.bed"

extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric",
                          qValue = "numeric", peak = "integer")

DAPseq_MACS2_narrowPeak <- import(file_narrowPeak, format = "BED",
                           extraCols = extraCols_narrowPeak)

DAPseq_MACS2_annotation <- annotatePeakInBatch(MACS2_narrowPeak,ensembl,"TSS",output = "nearestBiDirectionalPromoters",bindingRegion = c(-10000, 10000))

IE_DEGs<- as.matrix(read.table('~/Documents/Wheat/Transcriptomics_analysis/RKD_2021/data/3_aligned/DEGs/Immature embryo/Upregulated/IE_TM_vs_WE_upregulated_DEG_IDs_txt'))
ChIP_DAP_RKD_MACS2_gene_IDs <- DAPseq_MACS2_annotation$feature[RKD_ChIP_Wheat_MACS2_annotation$feature %in% DAPseq_MACS2_annotation$feature]
ChIP_DAP_RNAseq_RKD_MACS2_gene_IDs <- ChIP_DAP_RKD_MACS2_gene_IDs[ChIP_DAP_RKD_MACS2_gene_IDs %in% IE_DEGs]
ChIP_DAP_RNAseq_RKD_MACS2_gene_IDs

DAP_RNAseq_RKD_MACS2_gene_IDs <- DAP_MACS_10kb[DAP_MACS_10kb %in% IE_DEGs]

write.table(DAP_RNAseq_RKD_MACS2_gene_IDs,'DAPseq_RNAseq_MACS2_peaks.txt',quote=F,sep='\t')
write.table(DAPseq_MACS2_annotation,"DAPseq_MACS2_peaks.txt",quote=F,sep="\t")
write.table(DAPseq_MACS2_annotation$feature,"DAPseq_MACS2_gene_IDs.txt",quote=F,row.names = F,col.names = F)
write.table(ChIP_DAP_RKD_MACS2_gene_IDs,"ChIP_DAPseq_MACS2_gene_IDs.txt",quote=F,row.names = F,col.names = F)
write.table(ChIP_DAP_RNAseq_RKD_MACS2_gene_IDs,"ChIP_DAP_RNAseq_MACS2_gene_IDs.txt",quote=F,row.names = F,col.names = F)

####################################
### Overlap between DAP and ChIP ###
####################################

DAPseq<-unique(DAPseq_MACS2_narrowPeak)
ChIPseq<-unique(ChIPseq_MACS2_narrowPeak)

Overlap_DAP_ChIP<- findOverlapsOfPeaks(DAPseq,ChIPseq)

makeVennDiagram(Overlap_DAP_ChIP,Txdb=txdb,by="feature")

####################################
### Enriched heatmap - ChIPseq ###
####################################

genes <- transcripts(txdb)

targets <- genes
ExtendSize <- 1000
genes(txdb)

BIGWIG_ChIP_1<-rtracklayer::import('~/Documents/ChIPseq/BIGWIG/Wh_Cut_RKD2_RKD2_1S8_R1_001.rmdup.sorted.all.bigwig',format="BigWig")

mat1<-normalizeToMatrix(signal=BIGWIG,
                        target=targets,
                        background=0,
                        keep= c(0,0.99),
                        mean_mode="w0",
                        value_column="score",
                        extend=ExtendSize)


normalise_to_matrix <- function(signal, target, extend = 2000, w = max(extend)/50,
                                meth = FALSE, bed = FALSE,
                                value_column = if(!bed) {"score"}, # trick to pass a NULL if bed is TRUE
                                mapping_column = NULL, background = NA,
                                mean_mode = case_when(meth ~ "absolute", bed  ~ "coverage", TRUE ~ "w0"),
                                include_target = any(width(target) > 1),
                                target_ratio = median(width(target))/(extend + median(width(target)) + extend),
                                k = min(c(20, min(width(target)))), smooth = FALSE, smooth_fun = slider_median,
                                keep = c(0, 0.99), flip_upstream = FALSE) {
  normalizeToMatrix(
    signal = signal, target = target, extend = extend, w = w, value_column = value_column, 
    mapping_column = mapping_column, background = background, mean_mode = mean_mode, 
    include_target = include_target, target_ratio = target_ratio, k = k, 
    smooth = smooth, smooth_fun = smooth_fun, keep = keep, flip_upstream = flip_upstream
  )
}

head<-head(BIGWIG_ChIP_1,n=100)
mat1<-normalise_to_matrix(head, targets)

mat1<-normalizeToMatrix(BIGWIG_ChIP_1,targets,extend=10000)
mat2<-normalizeToMatrix(BIGWIG_DAP, targets, extend=10000)

mat1_means <- colMeans(mat1)
matrix<-cbind(as.matrix(mat1_means),c(1:ncol(mat1))*200)
write.table(mat_means)
colnames(matrix) <- c("Score","Position")

ggplot(as.data.frame(matrix),aes(x = Position,y = Score)) +
  geom_line() +
  scale_x_continuous(breaks = c(200,9200,10200,13800,14800,23600),labels = c("-10kb","-1kb", "TSS","TES","1kb","10kb")) +
  geom_vline(xintercept = c(9200,10200,13800,14800),linetype = "dashed") +
  labs(x = "Position",y = "Score") +
  theme_cowplot()

####################################
### Enriched heatmap - DAPseq ###
####################################

genes <- transcripts(txdb)

targets <- genes
ExtendSize <- 1000
targets.extended<- resize(targets, fix = "center", width = ExtendSize*2)
genes(txdb)

BIGWIG_DAP<-rtracklayer::import('~/Documents/DAPseq/Wheat/BIGWIG/TaRKD2_30.rmdup.chr.all.bw',format="BigWig")

normalise_to_matrix <- function(signal, target, extend = 2000, w = max(extend)/50,
                                meth = FALSE, bed = FALSE,
                                value_column = if(!bed) {"score"}, # trick to pass a NULL if bed is TRUE
                                mapping_column = NULL, background = NA,
                                mean_mode = case_when(meth ~ "absolute", bed  ~ "coverage", TRUE ~ "w0"),
                                include_target = any(width(target) > 1),
                                target_ratio = median(width(target))/(extend + median(width(target)) + extend),
                                k = min(c(20, min(width(target)))), smooth = FALSE, smooth_fun = slider_median,
                                keep = c(0, 0.99), flip_upstream = FALSE) {
  normalizeToMatrix(
    signal = signal, target = target, extend = extend, w = w, value_column = value_column, 
    mapping_column = mapping_column, background = background, mean_mode = mean_mode, 
    include_target = include_target, target_ratio = target_ratio, k = k, 
    smooth = smooth, smooth_fun = smooth_fun, keep = keep, flip_upstream = flip_upstream
  )
}

mat1<-normalise_to_matrix(BIGWIG_DAP, targets)

mat2_means <- colMeans(mat2)
matrix<-cbind(as.matrix(mat2_means),c(1:ncol(mat1))*200)

colnames(matrix) <- c("Score","Position")

ggplot(as.data.frame(matrix),aes(x = Position,y = Score)) +
  geom_line() +
  scale_x_continuous(breaks = c(200,9200,10200,13800,14800,23600),labels = c("-10kb","-1kb", "TSS","TES","1kb","10kb")) +
  geom_vline(xintercept = c(9200,10200,13800,14800),linetype = "dashed") +
  labs(x = "Position",y = "Score") +
  theme_cowplot()

####################################
### Expected Dist ###
####################################

# Expected distribution calculation

exons <- as.data.frame(exons(txdb, columns=NULL))
fiveUTRs <-  as.data.frame(fiveUTRsByTranscript(txdb))
threeUTRs <-  as.data.frame(threeUTRsByTranscript(txdb))
introns <-  as.data.frame(intronsByTranscript(txdb))
Downstream <-  as.data.frame(getPromoters(txdb,downstream=10000,upstream=0))
Upstream <-  as.data.frame(getPromoters(txdb,downstream=0,upstream=10000))

exon_size<-sum(exons$width)
fiveUTRs_size<-sum(fiveUTRs$width)
threeUTRs_size<-sum(threeUTRs$width)
introns_size<-sum(introns$width)
Downstream_size<-sum(Downstream$width)
Upstream_size<-sum(Upstream$width)
Intergenic_size<-genome_length-exon_size-fiveUTRs_size-threeUTRs_size-introns_size-Downstream_size-Upstream_size


exon_percentage<-(exon_size/genome_length)*100
intron_percentage<-(introns_size/genome_length)*100
fiveUTRs_percentage<-(fiveUTRs_size/genome_length)*100
threeUTRs_percentage<-(threeUTRs_size/genome_length)*100
Upstream_percentage<-(Upstream_size/genome_length)*100
Downstream_percentage<-(Downstream_size/genome_length)*100
Intergenic_percentage<-(Intergenic_size/genome_length)*100

Expected <- c(Upstream_percentage,Downstream_percentage,fiveUTRs_percentage,threeUTRs_percentage,
              exon_percentage,intron_percentage,Intergenic_percentage)

# Feature distribution calculation

library(ggplot2)


ChIPseq_aCR <- assignChromosomeRegion(ChIPseq_MACS2_narrowPeak,nucleotideLevel = F,
                              precedence = c("Promoters","immediateDownstream","fiveUTRs",
                                             "threeUTRs","Exons","Introns"),
                              immediate.downstream.cutoff=10000,
                              proximal.promoter.cutoff=10000,
                              TxDb=txdb)

DAPseq_aCR <- assignChromosomeRegion(DAPseq_MACS2_narrowPeak,nucleotideLevel = F,
                                     precedence = c("Promoters","immediateDownstream","fiveUTRs",
                                                    "threeUTRs","Exons","Introns"),
                                     immediate.downstream.cutoff=10000,
                                     proximal.promoter.cutoff=10000,
                                     TxDb=txdb)

ChIPseq <- as.vector(ChIPseq_aCR$percentage)
features<-c("Promoters","Downstream","fiveUTRs","threeUTRs","Exons","Introns","Intergenic")
ChIPseq <- as.data.frame(cbind(ChIPseq,features))

DAPseq <- as.vector(DAPseq_aCR$percentage)
features<-c("Promoters","Downstream","fiveUTRs","threeUTRs","Exons","Introns","Intergenic")
DAPseq <- as.data.frame(cbind(DAPseq,features))

colnames(DAPseq)<-c("percentage","feature")
colnames(ChIPseq)<-c("percentage","feature")
colnames(Expected_aCR)<-c("percentage","feature")

Data<-rbind(DAPseq,ChIPseq,Expected_aCR)
Data$group<-c(rep("DAPseq",7),rep("ChIPseq",7),rep("Expected",7))
Data$percentage <- as.numeric(Data$percentage)
library(cowplot)
ggplot(Data,aes(x=group,y=percentage,fill=feature))+geom_bar(position="stack",stat="identity")+theme_cowplot()


DAP.p.values<- as.numeric(DAPseq$percentage)/as.numeric(Expected_aCR$percentage)                               
ChIP.p.values <- as.numeric(ChIPseq$percentage)/as.numeric(Expected_aCR$percentage)

heatmap <- cbind(DAP.p.values,ChIP.p.values)
rownames(heatmap)<-DAPseq$feature

library(RColorBrewer)

library(pheatmap)
col.pal <- RColorBrewer::brewer.pal(9, "Reds")
heatmap<-heatmap[c(1:6),]
colnames(heatmap)<-c("DAPseq","ChIPseq")
pheatmap(log(heatmap), scale = "column",border_color=NA,show_rownames = T,cluster_cols = F, main = 'Observed vs Expected distribution',color=col.pal)

#No intergenic
Data <- Data[!Data$feature=="Intergenic",]

Data$group <- factor(Data$group,levels = c("Expected","DAPseq","ChIPseq"))

ggplot(Data,aes(x=group,y=percentage,fill=feature))+geom_bar(position="fill",stat="identity")+theme_cowplot()

#No intergenic No chIP

DAP.p.values<- as.numeric(DAPseq$percentage)/as.numeric(Expected_aCR$percentage)                               
ChIP.p.values <- as.numeric(ChIPseq$percentage)/as.numeric(Expected_aCR$percentage)

heatmap <- as.matrix(DAP.p.values)
rownames(heatmap)<-DAPseq$feature
heatmap<-heatmap[c(1:6),]

col.pal <- RColorBrewer::brewer.pal(9, "Reds")
pheatmap(log(heatmap), scale = "column",border_color=NA,show_rownames = T,cluster_cols = F, main = 'DAPseq vs Expected distribution',color=col.pal)

#No intergenic
Data <- Data[!Data$feature=="Intergenic",]
Data <- Data[!Data$group=="ChIPseq",]

ggplot(Data,aes(x=group,y=percentage,fill=feature))+geom_bar(position="fill",stat="identity")+theme_cowplot()

####################################
### ChIPseeker - ChIPseq ###
#################################### 

library(clusterProfiler)
covplot(ChIPseq_MACS2_narrowPeak, weightCol="score")

promoter <- getPromoters(TxDb=txdb, upstream=10000, downstream=10000)
tagMatrix <- getTagMatrix(ChIPseq_MACS2_narrowPeak, windows=promoter)
tagHeatmap(tagMatrix, xlim=c(-10000, 10000), color="red")

plotAvgProf(tagMatrix, xlim=c(-10000, 10000),
            xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")

peakAnno_ChIP <- annotatePeak(ChIPseq_MACS2_narrowPeak,tssRegion=c(-10000, 10000),
                               TxDb=txdb)
plotAnnoPie(peakAnno_ChIP)

plotAnnoBar(peakAnno_ChIP)

vennpie(peakAnno_ChIP)

upsetplot(peakAnno_ChIP)

plotDistToTSS(peakAnno_ChIP,
              title="Distribution of ChIPseq-binding loci\nrelative to TSS")

####################################
### ChIPseeker - DAPseq ###
#################################### 

library(clusterProfiler)
covplot(DAPseq_MACS2_narrowPeak, title = "DAP Peaks over Chromosomes", weightCol="score")

promoter <- getPromoters(TxDb=txdb, upstream=10000, downstream=10000)
tagMatrix <- getTagMatrix(DAPseq_MACS2_narrowPeak, windows=promoter)
tagHeatmap(tagMatrix, xlim=c(-10000, 10000), color="red")

plotAvgProf(tagMatrix, xlim=c(-10000, 10000),
            xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")

peakAnno_DAP <- annotatePeak(DAPseq_MACS2_narrowPeak,tssRegion=c(-10000, 10000),
                              TxDb=txdb)

plotAnnoPie(peakAnno_DAP)

plotAnnoBar(peakAnno_DAP)

vennpie(peakAnno_DAP)

upsetplot(peakAnno_DAP)

plotDistToTSS(peakAnno_DAP,
              title="Distribution of DAPseq-binding loci\nrelative to TSS")


####################################################################
library(ggvenn)
ChIPseq_10kb <- RKD_ChIP_Wheat_MACS2_annotation$feature
DAPseq_10kb <- DAPseq_MACS2_annotation$feature

x<-list(ChIPseq_10kb = ChIPseq_10kb,
        DAPseq_10kb = DAPseq_10kb)

ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF"),
  stroke_size = 0.5, set_name_size = 4
)

Upregulated_DEGs <- rownames(read.delim("~/Documents/Wheat/Transcriptomics_analysis/RKD_2021/data/3_aligned/DEGs/Immature embryo/Upregulated/TM_vs_NE_upregulated_DEGs.txt",sep=' '))
Downregulated_DEGs <- rownames(read.delim("~/Documents/Wheat/Transcriptomics_analysis/RKD_2021/data/3_aligned/DEGs/Immature embryo/Downregulated/IE_TM_vs_WE_downregulated_DEGs.txt",sep=' '))

ChIPseq_DAPseq_10kb <- ChIPseq_10kb[ChIPseq_10kb %in% DAPseq_10kb]

y <- list(DAPseq_10kb = DAPseq_10kb,
          Upregulated_DEGs = Upregulated_DEGs,
          Downregulated_DEGs = Downregulated_DEGs)

DAPseq_Up_10kb <- DAPseq_10kb[DAPseq_10kb %in% Upregulated_DEGs]
DAPseq_Down_10kb <- DAPseq_10kb[DAPseq_10kb %in% Downregulated_DEGs]

write.table(DAPseq_Up_10kb,"ChIPseq_DAPseq_Up_10kb_IDs.txt",quote=F,row.names = F,col.names = F)
write.table(DAPseq_Down_10kb,"ChIPseq_DAPseq_Down_10kb_IDs.txt",quote=F,row.names = F,col.names = F)


ggvenn(
  y, 
  fill_color = c("#0073C2FF", "#EFC000FF","#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
)

grid.newpage()
grid.draw(vp)

vp = venn.diagram(list(ChIP_MACS_10kb = ChIP_MACS_10kb,
                       DAP_MACS_10kb = DAP_MACS_10kb),
                  fill = c("red", "blue"),
                  alpha = c(0.5, 0.5), cex = 2, lty =2, 
                  filename = NULL)

grid.newpage()
grid.draw(vp)