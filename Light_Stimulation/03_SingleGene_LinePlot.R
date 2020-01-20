suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(ggforce))
suppressPackageStartupMessages(library(reshape2))

load("output/GeneSets_LightSim_Deg.RData")
RPKM=read.table("input/RPKM_LighSim_EXON.txt",sep="\t",header=T)
ord=c("CT17_1","CT17_1s","CT17_2","CT17_2s","Min30_1","Min30_1s","Min30_2","Min30_2s","Hour1_1","Hour1_2","Hour3_1","Hour3_2","Hour6_1","Hour6_2","CT23_1","CT23_2")
RPKM=RPKM[,match(ord,names(RPKM))]
cols <- c("0","0","0","0","0.5","0.5","0.5","0.5","1","1","3","3","6","6","CT23_1","CT23_2")
colnames(RPKM)=cols
RPKM$CT23_1 <- NULL
RPKM$CT23_2 <- NULL
RPKM$Rows <- rownames(RPKM)

df <- melt(RPKM)
df$variable <- as.numeric(as.character(df$variable))
df$log <- log2(df$value + 1)
df$Rows <- as.factor(df$Rows)
df$Labels <- df$Rows
df <- droplevels(df)

# CT17 vs 30 min
up <- GeneSets_LightSim_Deg[[1]][GeneSets_LightSim_Deg[[1]]$log2FoldChange < 0,]
down <- GeneSets_LightSim_Deg[[1]][GeneSets_LightSim_Deg[[1]]$log2FoldChange > 0,]
up <- up[order(up$log2FoldChange),]
down <- down[order(-down$log2FoldChange),]


filtUp <- df[df$Rows %in% rownames(up),]
filtUp <- droplevels(filtUp)
filtUp$Rows=factor(filtUp$Rows, levels = rownames(up))

filtDown <- df[df$Rows %in% rownames(down),]
filtDown <- droplevels(filtDown)
filtDown$Rows=factor(filtDown$Rows, levels = rownames(down))

ii <- c(0,0.5,1,2,3,4,5,6)
n_pages <- ceiling(length(levels(filtUp$Rows))/ 9)
pdf("visualizations/LinePlot_CT17vs30min_INDUCED_logRPKM.pdf",width=6,height=6,useDingbats=FALSE)
for (i in seq_len(n_pages)) {
print(ggplot(filtUp, aes(variable, log)) + 
geom_line(color = "steelblue", size=0.5,linetype = "dashed") +
geom_point(size=1)+
theme_classic()+
labs(title = "30 min: Induced")+
xlab("")+ 
ylab("log2(RPKM+1)") +
theme(legend.position="none")+
scale_x_continuous(breaks=ii) + 
theme(axis.text.x = element_text(angle = 45, hjust = 1,size=8),axis.text.y = element_text(size=8))+
facet_wrap_paginate(~Rows, ncol = 3, nrow = 3, scales = "free", page = i))
}
dev.off()

n_pages <- ceiling(length(levels(filtDown$Rows))/ 9)
pdf("visualizations/LinePlot_CT17vs30min_REDUCED_logRPKM.pdf",width=6,height=6,useDingbats=FALSE)
for (i in seq_len(n_pages)) {
print(ggplot(filtDown, aes(variable, log)) + 
geom_line(color = "steelblue", size=0.5,linetype = "dashed") +
geom_point(size=1)+
theme_classic()+
labs(title = "30 min: Reduced")+
xlab("")+ 
ylab("log2(RPKM+1)") +
theme(legend.position="none")+
scale_x_continuous(breaks=ii) + 
theme(axis.text.x = element_text(angle = 45, hjust = 1,size=8),axis.text.y = element_text(size=8))+
facet_wrap_paginate(~Rows, ncol = 3, nrow = 3, scales = "free", page = i))
}
dev.off()


# CT17 vs 1 hr
up <- GeneSets_LightSim_Deg[[2]][GeneSets_LightSim_Deg[[2]]$log2FoldChange < 0,]
down <- GeneSets_LightSim_Deg[[2]][GeneSets_LightSim_Deg[[2]]$log2FoldChange > 0,]
up <- up[order(up$log2FoldChange),]
down <- down[order(-down$log2FoldChange),]


filtUp <- df[df$Rows %in% rownames(up),]
filtUp <- droplevels(filtUp)
filtUp$Rows=factor(filtUp$Rows, levels = rownames(up))

filtDown <- df[df$Rows %in% rownames(down),]
filtDown <- droplevels(filtDown)
filtDown$Rows=factor(filtDown$Rows, levels = rownames(down))

n_pages <- ceiling(length(levels(filtUp$Rows))/ 12)
pdf("visualizations/LinePlot_CT17vs1hr_INDUCED_logRPKM.pdf",width=6,height=6,useDingbats=FALSE)
for (i in seq_len(n_pages)) {
print(ggplot(filtUp, aes(variable, log)) + 
geom_line(color = "steelblue", size=0.5,linetype = "dashed") +
geom_point(size=1)+
theme_classic()+
labs(title = "1 hr: Induced")+
xlab("")+ 
ylab("log2(RPKM+1)") +
theme(legend.position="none")+
scale_x_continuous(breaks=ii) + 
theme(axis.text.x = element_text(angle = 45, hjust = 1,size=8),axis.text.y = element_text(size=8))+
facet_wrap_paginate(~Rows, ncol = 3, nrow = 3, scales = "free", page = i))
}
dev.off()

n_pages <- ceiling(length(levels(filtDown$Rows))/ 12)
pdf("visualizations/LinePlot_CT17vs1hr_REDUCED_logRPKM.pdf",width=6,height=6,useDingbats=FALSE)
for (i in seq_len(n_pages)) {
print(ggplot(filtDown, aes(variable, log)) + 
geom_line(color = "steelblue", size=0.5,linetype = "dashed") +
geom_point(size=1)+
theme_classic()+
labs(title = "1 hr: Reduced")+
xlab("")+ 
ylab("log2(RPKM+1)") +
theme(legend.position="none")+
scale_x_continuous(breaks=ii) + 
theme(axis.text.x = element_text(angle = 45, hjust = 1,size=8),axis.text.y = element_text(size=8))+
facet_wrap_paginate(~Rows, ncol = 3, nrow = 3, scales = "free", page = i))
}
dev.off()

# CT17 vs 3 hr
up <- GeneSets_LightSim_Deg[[3]][GeneSets_LightSim_Deg[[3]]$log2FoldChange < 0,]
down <- GeneSets_LightSim_Deg[[3]][GeneSets_LightSim_Deg[[3]]$log2FoldChange > 0,]
up <- up[order(up$log2FoldChange),]
down <- down[order(-down$log2FoldChange),]


filtUp <- df[df$Rows %in% rownames(up),]
filtUp <- droplevels(filtUp)
filtUp$Rows=factor(filtUp$Rows, levels = rownames(up))

filtDown <- df[df$Rows %in% rownames(down),]
filtDown <- droplevels(filtDown)
filtDown$Rows=factor(filtDown$Rows, levels = rownames(down))

n_pages <- ceiling(length(levels(filtUp$Rows))/ 12)
pdf("visualizations/LinePlot_CT17vs3hr_INDUCED_logRPKM.pdf",width=6,height=6,useDingbats=FALSE)
for (i in seq_len(n_pages)) {
print(ggplot(filtUp, aes(variable, log)) + 
geom_line(color = "steelblue", size=0.5,linetype = "dashed") +
geom_point(size=1)+
theme_classic()+
labs(title = "3 hr: Induced")+
xlab("")+ 
ylab("log2(RPKM+1)") +
theme(legend.position="none")+
scale_x_continuous(breaks=ii) + 
theme(axis.text.x = element_text(angle = 45, hjust = 1,size=8),axis.text.y = element_text(size=8))+
facet_wrap_paginate(~Rows, ncol = 3, nrow = 3, scales = "free", page = i))
}
dev.off()

n_pages <- ceiling(length(levels(filtDown$Rows))/ 12)
pdf("visualizations/LinePlot_CT17vs3hr_REDUCED_logRPKM.pdf",width=6,height=6,useDingbats=FALSE)
for (i in seq_len(n_pages)) {
print(ggplot(filtDown, aes(variable, log)) + 
geom_line(color = "steelblue", size=0.5,linetype = "dashed") +
geom_point(size=1)+
theme_classic()+
labs(title = "3 hr: Reduced")+
xlab("")+ 
ylab("log2(RPKM+1)") +
theme(legend.position="none")+
scale_x_continuous(breaks=ii) + 
theme(axis.text.x = element_text(angle = 45, hjust = 1,size=8),axis.text.y = element_text(size=8))+
facet_wrap_paginate(~Rows, ncol = 3, nrow = 3, scales = "free", page = i))
}
dev.off()


# CT17 vs 6 hr
up <- GeneSets_LightSim_Deg[[4]][GeneSets_LightSim_Deg[[4]]$log2FoldChange < 0,]
down <- GeneSets_LightSim_Deg[[4]][GeneSets_LightSim_Deg[[4]]$log2FoldChange > 0,]

up <- up[rownames(up)!= rownames(GeneSets_LightSim_Deg[[5]]),]
down <- down[rownames(down)!= rownames(GeneSets_LightSim_Deg[[5]]),]

up <- up[order(up$log2FoldChange),]
down <- down[order(-down$log2FoldChange),]


filtUp <- df[df$Rows %in% rownames(up),]
filtUp <- droplevels(filtUp)
filtUp$Rows=factor(filtUp$Rows, levels = rownames(up))

filtDown <- df[df$Rows %in% rownames(down),]
filtDown <- droplevels(filtDown)
filtDown$Rows=factor(filtDown$Rows, levels = rownames(down))

n_pages <- ceiling(length(levels(filtUp$Rows))/ 12)
pdf("visualizations/LinePlot_CT17vs6hr_INDUCED_logRPKM.pdf",width=6,height=6,useDingbats=FALSE)
for (i in seq_len(n_pages)) {
print(ggplot(filtUp, aes(variable, log)) + 
geom_line(color = "steelblue", size=0.5,linetype = "dashed") +
geom_point(size=1)+
theme_classic()+
labs(title = "6 hr: Induced")+
xlab("")+ 
ylab("log2(RPKM+1)") +
theme(legend.position="none")+
scale_x_continuous(breaks=ii) + 
theme(axis.text.x = element_text(angle = 45, hjust = 1,size=8),axis.text.y = element_text(size=8))+
facet_wrap_paginate(~Rows, ncol = 3, nrow = 3, scales = "free", page = i))
}
dev.off()

n_pages <- ceiling(length(levels(filtDown$Rows))/ 12)
pdf("visualizations/LinePlot_CT17vs6hr_REDUCED_logRPKM.pdf",width=6,height=6,useDingbats=FALSE)
for (i in seq_len(n_pages)) {
print(ggplot(filtDown, aes(variable, log)) + 
geom_line(color = "steelblue", size=0.5,linetype = "dashed") +
geom_point(size=1)+
theme_classic()+
labs(title = "6 hr: Reduced")+
xlab("")+ 
ylab("log2(RPKM+1)") +
theme(legend.position="none")+
scale_x_continuous(breaks=ii) + 
theme(axis.text.x = element_text(angle = 45, hjust = 1,size=8),axis.text.y = element_text(size=8))+
facet_wrap_paginate(~Rows, ncol = 3, nrow = 3, scales = "free", page = i))
}
dev.off()

