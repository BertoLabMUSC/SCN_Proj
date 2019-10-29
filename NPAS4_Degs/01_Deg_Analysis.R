# Library to load for this step
library(DESeq2)
library(xlsx)

rm(list=ls())
RPKM=read.table("input/RPKM_NPAS4_EXON.txt",sep="\t",header=T)
COUNT=read.table("input/COUNT_NPAS4_EXON.txt",sep="\t",header=T)

#WDvsWL
WDvsWL=COUNT[c(7,8,9,10,11,12)]
tmp <- RPKM[c(7,8,9,10,11,12)]
first=apply(tmp, 1, function(x) (all(x[1:3] >= 0.5) | all(x[4:6] >= 0.5)))
WDvsWL <- WDvsWL[first,]

DESIGN=data.frame(row.names = colnames(WDvsWL), Treatment=c(rep("DD",3),rep("LI",3)))
dds <- DESeqDataSetFromMatrix(countData = WDvsWL,colData = DESIGN,design = ~ Treatment)
dds$Treatment <- factor(dds$Treatment, levels=c("DD", "LI"))
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds, full=design(dds), betaPrior=TRUE)
res <- as.data.frame(results(dds,contrast=c("Treatment","DD","LI"),cooksCutoff=FALSE))
res$log2FoldChange <- -1*res$log2FoldChange
res$absLog=abs(res$log2FoldChange)
write.table(res,"OUTPUTS_woSVA/DESeq_NPAS4_WDvsWL.txt",sep="\t",quote=F)

#KDvsKL
rm(list=ls())
RPKM=read.table("input/RPKM_NPAS4_EXON.txt",sep="\t",header=T)
COUNT=read.table("input/COUNT_NPAS4_EXON.txt",sep="\t",header=T)
KDvsKL=COUNT[c(1,2,3,4,5,6)]
tmp <- RPKM[c(1,2,3,4,5,6)]
first=apply(tmp, 1, function(x) (all(x[1:3] >= 0.5) | all(x[4:6] >= 0.5)))
KDvsKL <- KDvsKL[first,]

DESIGN=data.frame(row.names = colnames(KDvsKL), Treatment=c(rep("DD",3),rep("LI",3)))
dds <- DESeqDataSetFromMatrix(countData = KDvsKL,colData = DESIGN,design = ~ Treatment)
dds$Treatment <- factor(dds$Treatment, levels=c("DD", "LI"))
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds, full=design(dds), betaPrior=TRUE)
res <- as.data.frame(results(dds,contrast=c("Treatment","DD","LI"),cooksCutoff=FALSE))
res$log2FoldChange <- -1*res$log2FoldChange
res$absLog=abs(res$log2FoldChange)
write.table(res,"output/DESeq_NPAS4_KDvsKL.txt",sep="\t",quote=F)

#WDvsKD
rm(list=ls())
RPKM=read.table("input/RPKM_NPAS4_EXON.txt",sep="\t",header=T)
COUNT=read.table("input/COUNT_NPAS4_EXON.txt",sep="\t",header=T)
WDvsKD=COUNT[c(7,8,9,1,2,3)]
tmp <- RPKM[c(7,8,9,1,2,3)]
first=apply(tmp, 1, function(x) (all(x[1:3] >= 0.5) | all(x[4:6] >= 0.5)))
WDvsKD <- WDvsKD[first,]

DESIGN=data.frame(row.names = colnames(WDvsKD), Treatment=c(rep("DD",3),rep("LI",3)))
dds <- DESeqDataSetFromMatrix(countData = WDvsKD,colData = DESIGN,design = ~ Treatment)
dds$Treatment <- factor(dds$Treatment, levels=c("DD", "LI"))
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds, full=design(dds), betaPrior=TRUE)
res <- as.data.frame(results(dds,contrast=c("Treatment","DD","LI"),cooksCutoff=FALSE))
res$log2FoldChange <- -1*res$log2FoldChange
res$absLog=abs(res$log2FoldChange)
write.table(res,"output/DESeq_NPAS4_WDvsKD.txt",sep="\t",quote=F)


#WLvsKL
rm(list=ls())
RPKM=read.table("input/RPKM_NPAS4_EXON.txt",sep="\t",header=T)
COUNT=read.table("input/COUNT_NPAS4_EXON.txt",sep="\t",header=T)
WLvsKL=COUNT[c(10,11,12,4,5,6)]
tmp <- RPKM[c(10,11,12,4,5,6)]
first=apply(tmp, 1, function(x) (all(x[1:3] >= 0.5) | all(x[4:6] >= 0.5)))
WLvsKL <- WLvsKL[first,]

DESIGN=data.frame(row.names = colnames(WLvsKL), Treatment=c(rep("DD",3),rep("LI",3)))
dds <- DESeqDataSetFromMatrix(countData = WLvsKL,colData = DESIGN,design = ~ Treatment)
dds$Treatment <- factor(dds$Treatment, levels=c("DD", "LI"))
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds, full=design(dds), betaPrior=TRUE)
res <- as.data.frame(results(dds,contrast=c("Treatment","DD","LI"),cooksCutoff=FALSE))
res$log2FoldChange <- -1*res$log2FoldChange
res$absLog=abs(res$log2FoldChange)
write.table(res,"output/DESeq_NPAS4_WLvsKL.txt",sep="\t",quote=F)


## 
KDvsKL=read.table("output/DESeq_NPAS4_KDvsKL.txt")
WDvsKD=read.table("output/DESeq_NPAS4_WDvsKD.txt")
WDvsWL=read.table("output/DESeq_NPAS4_WDvsWL.txt")
WLvsKL=read.table("output/DESeq_NPAS4_WLvsKL.txt")

KDvsKL=KDvsKL[KDvsKL$padj < 0.05,]
WDvsKD=WDvsKD[WDvsKD$padj < 0.05,]
WDvsWL=WDvsWL[WDvsWL$padj < 0.05,]
WLvsKL=WLvsKL[WLvsKL$padj < 0.05,]

write.xlsx(KDvsKL, file="output/NPAS4_DEGs.xlsx", sheetName="KDvsKL", row.names=TRUE)
write.xlsx(WDvsKD, file="output/NPAS4_DEGs.xlsx", sheetName="WDvsKD", append=TRUE, row.names=TRUE)
write.xlsx(WDvsWL, file="output/NPAS4_DEGs.xlsx", sheetName="WDvsWL", append=TRUE, row.names=TRUE)
write.xlsx(WLvsKL, file="output/NPAS4_DEGs.xlsx", sheetName="WLvsKL", append=TRUE, row.names=TRUE)
            
sessionInfo()
