# Libraries
suppressPackageStartupMessages(library(ggpubr))

df1 <- read.table("output/DESeq_NPAS4_WDvsWL.txt")
df1 <- df1[c(2,6)]
colnames(df1) <- c("logFC_WDWL","FDR_WDWL")

df2 <- read.table("output/DESeq_NPAS4_KDvsKL.txt")
df2 <- df2[c(2,6)]
colnames(df2) <- c("logFC_KDKL","FDR_KDKL")

x <- merge(df1,df2,by="row.names",all=FALSE)
sign <- x[x$FDR_WDWL < 0.05 | x$FDR_KDKL < 0.05,]
sign$Class <- ifelse(sign$FDR_WDWL < 0.05 & sign$FDR_KDKL < 0.05,"BOTH",
				ifelse(sign$FDR_WDWL < 0.05 & sign$FDR_KDKL > 0.05,"WDWLsp",
				 ifelse(sign$FDR_WDWL > 0.05 & sign$FDR_KDKL < 0.05,"KDKLsp","NA")))
sign2 <- sign[abs(sign$logFC_WDWL) > 0.3 | abs(sign$logFC_KDKL) > 0.3,]

pdf("FC_Comparison_WLWD_KLKD.pdf",width=5,height=5,useDingbats=FALSE)
ggscatter(sign2, x = "logFC_WDWL", y = "logFC_KDKL",
   color = "Class",size = 2,rug = FALSE,shape=1)+
geom_vline(xintercept = 0, colour = "grey60",linetype="dotted",size=1,alpha=0.5) + 
geom_hline(yintercept = 0, colour = "grey60",linetype="dotted",size=1,alpha=0.5) +
geom_abline(intercept = 0, colour = "grey60",linetype="dotted",size=1,alpha=0.5)+
xlim(-1.5,5.5)+
ylim(-1.5,5.5)+
theme(legend.position=c(0.15,0.8))+
ggtitle("Fold Change Comparison")+
xlab("log2FC WL/WD")+ 
ylab("log2FC KL/KD")
dev.off()







