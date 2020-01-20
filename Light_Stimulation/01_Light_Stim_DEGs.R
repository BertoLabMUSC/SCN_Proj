suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(sva))

dir.create("output/")

# Deg By Time
# CT17 vs 30 min
rm(list=ls())
RPKM=read.table("inputs/RPKM_LighSim_EXON.txt",sep="\t",header=T)
COUNT=read.table("inputs/COUNT_LighSim_EXON.txt",sep="\t",header=T)
ord=c("CT17_1","CT17_1s","CT17_2","CT17_2s","Min30_1","Min30_1s","Min30_2","Min30_2s","Hour1_1","Hour1_2","Hour3_1","Hour3_2","Hour6_1","Hour6_2","CT23_1","CT23_2")
RPKM=RPKM[,match(ord,names(RPKM))]
RPKM=RPKM[,1:14]

first=RPKM[c(1,2,3,4,5,6,7,8)]
first=apply(first, 1, function(x) (all(x[1:4] >= 0.5) | all(x[5:8] >= 0.5)))

COUNT=COUNT[first,match(ord,names(COUNT))]
tab=COUNT[c(1,2,3,4,5,6,7,8)]

DESIGN=data.frame(row.names = colnames(tab), Treatment=c(rep("DD",4),rep("LI",4)), Batch = as.factor(c(1,2,1,2,1,2,1,2)))
dds <- DESeqDataSetFromMatrix(countData = tab,colData = DESIGN,design = ~ Treatment)
dds$Treatment <- factor(dds$Treatment, levels=c("DD", "LI"))
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds, full=design(dds), betaPrior=TRUE)
res <- as.data.frame(results(dds,contrast=c("Treatment","DD","LI"),cooksCutoff=TRUE,independentFiltering = FALSE))
write.table(res,"output/DESeq_LIGHTSIM_CT17vs30min.txt",sep="\t",quote=F)

# Deg By Time
# CT17 vs 1hr
rm(list=ls())
RPKM=read.table("inputs/RPKM_LighSim_EXON.txt",sep="\t",header=T)
COUNT=read.table("inputs/COUNT_LighSim_EXON.txt",sep="\t",header=T)
ord=c("CT17_1","CT17_1s","CT17_2","CT17_2s","Min30_1","Min30_1s","Min30_2","Min30_2s","Hour1_1","Hour1_2","Hour3_1","Hour3_2","Hour6_1","Hour6_2","CT23_1","CT23_2")
RPKM=RPKM[,match(ord,names(RPKM))]
RPKM=RPKM[,1:14]

first=RPKM[c(1,2,3,4,9,10)]
first=apply(first, 1, function(x) (all(x[1:4] >= 0.5) | all(x[5:6] >= 0.5)))

COUNT=COUNT[first,match(ord,names(COUNT))]
tab=COUNT[c(1,2,3,4,9,10)]

DESIGN=data.frame(row.names = colnames(tab), Treatment=c(rep("DD",4),rep("LI",2)))
dds <- DESeqDataSetFromMatrix(countData = tab,colData = DESIGN,design = ~ Treatment)
dds$Treatment <- factor(dds$Treatment, levels=c("DD", "LI"))
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds, full=design(dds), betaPrior=TRUE)
res <- as.data.frame(results(dds,contrast=c("Treatment","DD","LI"),cooksCutoff=TRUE,independentFiltering = FALSE))
write.table(res,"output/DESeq_LIGHTSIM_CT17vs1hr.txt",sep="\t",quote=F)

# Deg By Time
# CT17 vs 3hr
rm(list=ls())
RPKM=read.table("inputs/RPKM_LighSim_EXON.txt",sep="\t",header=T)
COUNT=read.table("inputs/COUNT_LighSim_EXON.txt",sep="\t",header=T)
ord=c("CT17_1","CT17_1s","CT17_2","CT17_2s","Min30_1","Min30_1s","Min30_2","Min30_2s","Hour1_1","Hour1_2","Hour3_1","Hour3_2","Hour6_1","Hour6_2","CT23_1","CT23_2")
RPKM=RPKM[,match(ord,names(RPKM))]
RPKM=RPKM[,1:14]

first=RPKM[c(1,2,3,4,11,12)]
first=apply(first, 1, function(x) (all(x[1:4] >= 0.5) | all(x[5:6] >= 0.5)))

COUNT=COUNT[first,match(ord,names(COUNT))]
tab=COUNT[c(1,2,3,4,11,12)]

DESIGN=data.frame(row.names = colnames(tab), Treatment=c(rep("DD",4),rep("LI",2)))
dds <- DESeqDataSetFromMatrix(countData = tab,colData = DESIGN,design = ~ Treatment)
dds$Treatment <- factor(dds$Treatment, levels=c("DD", "LI"))
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds, full=design(dds), betaPrior=TRUE)
res <- as.data.frame(results(dds,contrast=c("Treatment","DD","LI"),cooksCutoff=TRUE,independentFiltering = FALSE))
write.table(res,"output/DESeq_LIGHTSIM_CT17vs3hr.txt",sep="\t",quote=F)


# Deg By Time
# CT17 vs 3hr
rm(list=ls())
RPKM=read.table("inputs/RPKM_LighSim_EXON.txt",sep="\t",header=T)
COUNT=read.table("inputs/COUNT_LighSim_EXON.txt",sep="\t",header=T)
ord=c("CT17_1","CT17_1s","CT17_2","CT17_2s","Min30_1","Min30_1s","Min30_2","Min30_2s","Hour1_1","Hour1_2","Hour3_1","Hour3_2","Hour6_1","Hour6_2","CT23_1","CT23_2")
RPKM=RPKM[,match(ord,names(RPKM))]
RPKM=RPKM[,1:14]

first=RPKM[c(1,2,3,4,13,14)]
first=apply(first, 1, function(x) (all(x[1:4] >= 0.5) | all(x[5:6] >= 0.5)))

COUNT=COUNT[first,match(ord,names(COUNT))]
tab=COUNT[c(1,2,3,4,13,14)]

DESIGN=data.frame(row.names = colnames(tab), Treatment=c(rep("DD",4),rep("LI",2)))
dds <- DESeqDataSetFromMatrix(countData = tab,colData = DESIGN,design = ~ Treatment)
dds$Treatment <- factor(dds$Treatment, levels=c("DD", "LI"))
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds, full=design(dds), betaPrior=TRUE)
res <- as.data.frame(results(dds,contrast=c("Treatment","DD","LI"),cooksCutoff=TRUE,independentFiltering = FALSE))
write.table(res,"output/DESeq_LIGHTSIM_CT17vs6hr.txt",sep="\t",quote=F)


# Deg By Time
# CT17 vs 3hr
rm(list=ls())
RPKM=read.table("inputs/RPKM_LighSim_EXON.txt",sep="\t",header=T)
COUNT=read.table("inputs/COUNT_LighSim_EXON.txt",sep="\t",header=T)
ord=c("CT17_1","CT17_1s","CT17_2","CT17_2s","Min30_1","Min30_1s","Min30_2","Min30_2s","Hour1_1","Hour1_2","Hour3_1","Hour3_2","Hour6_1","Hour6_2","CT23_1","CT23_2")
RPKM=RPKM[,match(ord,names(RPKM))]

first=RPKM[c(1,2,3,4,15,16)]
first=apply(first, 1, function(x) (all(x[1:4] >= 0.5) | all(x[5:6] >= 0.5)))

COUNT=COUNT[first,match(ord,names(COUNT))]
tab=COUNT[c(1,2,3,4,15,16)]

DESIGN=data.frame(row.names = colnames(tab), Treatment=c(rep("DD",4),rep("LI",2)))
dds <- DESeqDataSetFromMatrix(countData = tab,colData = DESIGN,design = ~ Treatment)
dds$Treatment <- factor(dds$Treatment, levels=c("DD", "LI"))
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds, full=design(dds), betaPrior=TRUE)
res <- as.data.frame(results(dds,contrast=c("Treatment","DD","LI"),cooksCutoff=TRUE,independentFiltering = FALSE))
write.table(res,"output/DESeq_LIGHTSIM_CT17vsCT23.txt",sep="\t",quote=F)


# Create Tables 
library(xlsx)
a=read.table("output/DESeq_LIGHTSIM_CT17vs30min.txt")
b=read.table("output/DESeq_LIGHTSIM_CT17vs1hr.txt")
c=read.table("output/DESeq_LIGHTSIM_CT17vs3hr.txt")
d=read.table("output/DESeq_LIGHTSIM_CT17vs6hr.txt")
e=read.table("output/DESeq_LIGHTSIM_CT17vsCT23.txt")

a=a[a$padj < 0.05 & abs(a$log2FoldChange) > 0.3,]
b=b[b$padj < 0.05 & abs(b$log2FoldChange) > 0.3,]
c=c[c$padj < 0.05 & abs(c$log2FoldChange) > 0.3,]
d=d[d$padj < 0.05 & abs(d$log2FoldChange) > 0.3,]
e=e[e$padj < 0.05 & abs(e$log2FoldChange) > 0.3,]


write.xlsx(a, file="output/table_S1_LightStim_DGE.xsls", sheetName="CT17vs30min", row.names=TRUE)
write.xlsx(b, file="output/table_S1_LightStim_DGE.xlsx", sheetName="CT17vs1hr", append=TRUE, row.names=TRUE)
write.xlsx(c, file="output/table_S1_LightStim_DGE.xlsx", sheetName="CT17vs3hr", append=TRUE, row.names=TRUE)
write.xlsx(d, file="output/table_S1_LightStim_DGE.xlsx", sheetName="CT17vs6hr", append=TRUE, row.names=TRUE)
write.xlsx(e, file="output/table_S1_LightStim_DGE.xlsx", sheetName="CT17vsCT23", append=TRUE, row.names=TRUE)

# Create Gene Set data 
CT17vs30min=read.table("output/DESeq_LIGHTSIM_CT17vs30min.txt")
CT17vs30min=CT17vs30min[CT17vs30min$padj < 0.05 & abs(CT17vs30min$log2FoldChange)>0.3,]

CT17vs1hr=read.table("output/DESeq_LIGHTSIM_CT17vs1hr.txt")
CT17vs1hr=CT17vs1hr[CT17vs1hr$padj < 0.05 & abs(CT17vs1hr$log2FoldChange)>0.3,]

CT17vs3hr=read.table("output/DESeq_LIGHTSIM_CT17vs3hr.txt")
CT17vs3hr=CT17vs3hr[CT17vs3hr$padj < 0.05 & abs(CT17vs3hr$log2FoldChange)>0.3,]

CT17vs6hr=read.table("DESeq_LIGHTSIM_CT17vs6hr.txt")
CT17vs6hr=CT17vs6hr[CT17vs6hr$padj < 0.05 & abs(CT17vs6hr$log2FoldChange)>0.3,]


CT17vsCT23=read.table("output/DESeq_LIGHTSIM_CT17vsCT23.txt")
CT17vsCT23=CT17vsCT23[CT17vsCT23$padj < 0.05 & abs(CT17vsCT23$log2FoldChange)>0.3,]

GeneSets_LightSim_Deg=list(CT17vs30min,CT17vs1hr,CT17vs3hr,CT17vs6hr,CT17vsCT23)
names(GeneSets_LightSim_Deg)=c("CT17vs30min","CT17vs1hr","CT17vs3hr","CT17vs6hr","CT17vsCT23")
save(GeneSets_LightSim_Deg,file="output/GeneSets_LightSim_Deg.RData")






