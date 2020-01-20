suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(ggforce))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(cowplot))

load("ouput/GeneSets_LightSim_Deg.RData")
RPKM=read.table("input/RPKM_LighSim_EXON.txt",sep="\t",header=T)
ord=c("CT17_1","CT17_1s","CT17_2","CT17_2s","Min30_1","Min30_1s","Min30_2","Min30_2s","Hour1_1","Hour1_2","Hour3_1","Hour3_2","Hour6_1","Hour6_2","CT23_1","CT23_2")
RPKM <- RPKM[,match(ord,colnames(RPKM))]
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

cool <- c("Egr1","Egr2","Egr3","Egr4",)
filt <- df[df$Rows %in% cool,]


ii <- c(0,0.5,1,2,3,4,5,6)
pdf("visualizations/CoolGuys_RPKM.pdf",width=3,height=3,useDingbats=FALSE)
plot1a <- ggplot(filt, aes(variable, value, color=Labels)) + 
geom_line(size=1) +
theme_classic()+
xlab("")+ 
ylab("log2(RPKM+1)") +
theme(legend.position="none")+
scale_x_continuous(breaks=ii) + 
theme(axis.text.x = element_text(angle = 45, hjust = 1,size=8),axis.text.y = element_text(size=8))
#ylim(0,40)
print(plot1a)
dev.off()


ii <- c(0,0.5,1,2,3,4,5,6)
pdf("visualizations/CoolGuys_RPKM.pdf",width=3,height=3,useDingbats=FALSE)
plot1a <- ggplot(filt, aes(variable, value, color=Labels)) + 
geom_line(size=1) +
theme_classic()+
xlab("")+ 
ylab("RPKM") +
theme(legend.position="none")+
scale_x_continuous(breaks=ii) + 
theme(axis.text.x = element_text(angle = 45, hjust = 1,size=8),axis.text.y = element_text(size=8))
#ylim(0,40)
print(plot1a)
dev.off()

cool <- c("Egr1","Egr2","Egr3","Egr4","Nr4a1","Nr4a2","Nr4a3")
filt <- df[df$Rows %in% cool,]


ii <- c(0,0.5,1,2,3,4,5,6)
pdf("visualizations/CoolGuys_EGRs.pdf",width=4,height=3,useDingbats=FALSE)
ggplot(filt, aes(variable, log, color=Labels)) + 
geom_line(size=1) +
theme_classic()+
xlab("")+ 
ylab("log2(RPKM+1)") +
theme(legend.position="right")+
scale_x_continuous(breaks=ii) + 
theme(axis.text.x = element_text(angle = 45, hjust = 1,size=8),axis.text.y = element_text(size=8))
#ylim(0,40)
dev.off()

cool <- c("Fos","Fosb","Fosl2","Junb","Jun")
filt <- df[df$Rows %in% cool,]


ii <- c(0,0.5,1,2,3,4,5,6)
pdf("visualizations/CoolGuys_AP-1.pdf",width=4,height=3,useDingbats=FALSE)
ggplot(filt, aes(variable, log, color=Labels)) + 
geom_line(size=1) +
theme_classic()+
xlab("")+ 
ylab("log2(RPKM+1)") +
theme(legend.position="right")+
scale_x_continuous(breaks=ii) + 
theme(axis.text.x = element_text(angle = 45, hjust = 1,size=8),axis.text.y = element_text(size=8))
dev.off()

cool <- c("Sik1","Dusp1","Dusp4","Ppp1r15a","Plk2")
filt <- df[df$Rows %in% cool,]


ii <- c(0,0.5,1,2,3,4,5,6)
pdf("visualizations/CoolGuys_KinPhos.pdf",width=4,height=3,useDingbats=FALSE)
ggplot(filt, aes(variable, log, color=Labels)) + 
geom_line(size=1) +
theme_classic()+
xlab("")+ 
ylab("log2(RPKM+1)") +
theme(legend.position="right")+
scale_x_continuous(breaks=ii) + 
theme(axis.text.x = element_text(angle = 45, hjust = 1,size=8),axis.text.y = element_text(size=8))
dev.off()


cool <- c("Per1","Per2","Per3","Cry1","Cry2","Bhlhe40","Bhlhe41")
filt <- df[df$Rows %in% cool,]


ii <- c(0,0.5,1,2,3,4,5,6)
pdf("visualizations/CoolGuys_Circ.pdf",width=4,height=3,useDingbats=FALSE)
ggplot(filt, aes(variable, log, color=Labels)) + 
geom_line(size=1) +
theme_classic()+
xlab("")+ 
ylab("log2(RPKM+1)") +
theme(legend.position="right")+
scale_x_continuous(breaks=ii) + 
theme(axis.text.x = element_text(angle = 45, hjust = 1,size=8),axis.text.y = element_text(size=8))
dev.off()
