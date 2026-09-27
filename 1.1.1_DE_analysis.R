library("DESeq2")
library("RColorBrewer")
library("ggplot2")
library("matrixStats")

path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(path)

## DESeq analysis ----
txi <- read.delim("0_Input_raw_counts.tsv", row.names=1, check.names = F)
TPM <- read.delim("0_Input_TPM.tsv", row.names=1, check.names = F)

phenoData <- read.delim("0_Input_Metadata.tsv", check.names=F, row.names=1) %>% filter(ToFilter == "N")

txi <- txi[, row.names(phenoData)]
TPM <- TPM[, row.names(phenoData)]

dds <- DESeqDataSetFromMatrix(txi,
                                   colData = phenoData,
                                   design = ~ Class)

dds <- DESeq(dds, parallel = T)

## Comparison vs WT
DEresEM <- results(dds, contrast=c("Class", "E", "M"))
DEresSM <- results(dds, contrast=c("Class", "S", "M"))
DEresSE <- results(dds, contrast=c("Class", "E", "S"))

write.table(DEresEM, "2.2.1_DE_results_E_vs_M.tsv", sep="\t", quote=FALSE, row.names=T)
write.table(DEresSM, "2.2.1_DE_results_S_vs_M.tsv", sep="\t", quote=FALSE, row.names=T)
write.table(DEresSE, "2.2.1_DE_results_E_vs_S.tsv", sep="\t", quote=FALSE, row.names=T)

## 2.1.1 Quality check PCA ----
library("ggfortify")
TPM <- read.delim("0_Input_TPM.tsv", row.names=1, check.names = F)
phenoData <- read.delim("0_Input_Metadata.tsv", check.names=F, row.names=1) %>% filter(ToFilter == "N")

phenoData$Class <- factor(phenoData$Class, levels=c("S", "M", "E"))

TPM <- TPM[, row.names(phenoData)]

TPM <- TPM[rowSums(TPM) > 0, ]
cdatat <- data.frame(t(TPM), check.names = F)
cdata.pca<-prcomp(log(cdatat+0.1,2), center=T, scale=T)

dat <- data.frame(cbind(phenoData, cdatat), check.names=F)

library("ggforce")

autoplot(cdata.pca, data = dat, colour = "Class", size=2, alpha=0.75, label=F, label.vjust=2, label.size=2) + 
  theme_bw() + 
  theme(legend.position="bottom") +
  ggforce::geom_mark_ellipse(aes(fill = Class, color = Class), n = 100) +
  labs(col="Class")

## 2.1.2 Pairwise heatmap ----
TPM <- read.delim("0_Input_TPM.tsv", row.names=1, check.names = F)
phenoData <- read.delim("0_Input_Metadata.tsv", check.names=F) %>% filter(ToFilter == "N")

TPM <- TPM[, phenoData$ID]

TPM <- TPM[rowSums(TPM) > 0, ]

sampleDists <- dist(t(TPM))

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- phenoData$ID
colnames(sampleDistMatrix) <- phenoData$ID
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

# Compute median and average expression
TPM <- read.delim("0_Input_TPM.tsv", row.names=1, check.names = F)
phenoData <- read.delim("0_Input_Metadata.tsv", check.names=F, row.names=1)  %>% filter(ToFilter == "N")
cov_class <- levels(factor(phenoData[,"Class"]))

TPM <- TPM[, row.names(phenoData)]

medt <- matrix(NA, ncol=length(cov_class), nrow=nrow(TPM)) 
avgt <- matrix(NA, ncol=length(cov_class), nrow=nrow(TPM))

colnames(medt) <- cov_class
colnames(avgt) <- cov_class

for(i in 1:length(cov_class)){
  
  subp <- subset(phenoData, phenoData$Class == cov_class[i])
  subc <- TPM[, which(colnames(TPM) %in% rownames(subp))]
  
  medt[,i] <- rowMedians(as.matrix(subc))
  avgt[,i] <- rowMeans(subc)
  
  colnames(medt)[i] <- paste0("Median ", cov_class[i])
  colnames(avgt)[i] <- paste0("Average ", cov_class[i])
  
}

res <- cbind(medt, avgt)
row.names(res) <- rownames(TPM)

write.table(res, "2.1_Expression_Level_median_avg.tsv", sep="\t", quote=FALSE, row.names=T)