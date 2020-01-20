suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(ggforce))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(cowplot))

dir.create("visualizations/")

load("output/GeneSets_LightSim_Deg.RData")
RPKM=read.table("inputs/RPKM_LighSim_EXON.txt",sep="\t",header=T)
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

# CT17 vs 30 min
up30 <- GeneSets_LightSim_Deg[[1]][GeneSets_LightSim_Deg[[1]]$log2FoldChange < 0,]
down30 <- GeneSets_LightSim_Deg[[1]][GeneSets_LightSim_Deg[[1]]$log2FoldChange > 0,]

filtUp30 <- df[df$Rows %in% rownames(up30),]
filtUp30 <- droplevels(filtUp30)

filtDown30 <- df[df$Rows %in% rownames(down30),]
filtDown30 <- droplevels(filtDown30)

ii <- c(0,0.5,1,2,3,4,5,6)
pdf("visualizations/MultiLine_30_INDUCED.pdf",width=3,height=4,useDingbats=FALSE)
plot1a <- ggplot(filtUp30, aes(variable, log, color=Labels)) + 
geom_line(alpha=0.5) +
theme_classic()+
#scale_color_manual(values=rep("grey",length(filtUp30$Labels)))+
stat_summary(fun.y = "mean", geom = "line", color = "black", size = 1)+
labs(title = "30 min: Induced", subtitle = "Dynamic")+
xlab("")+ 
ylab("log2(RPKM+1)") +
theme(legend.position="none")+
scale_x_continuous(breaks=ii) + 
theme(axis.text.x = element_text(angle = 45, hjust = 1,size=8),axis.text.y = element_text(size=8)) +
ylim(0,10)
print(plot1a)
dev.off()


ii <- c(0,0.5,1,2,3,4,5,6)
pdf("visualizations/MultiLine_30_REDUCED.pdf",width=3,height=4,useDingbats=FALSE)
plot1b <- ggplot(filtDown30, aes(variable, log, color=Labels)) + 
geom_line(alpha=0.5) +
theme_classic()+
#scale_color_manual(values=rep("grey",length(filtUp30$Labels)))+
stat_summary(fun.y = "mean", geom = "line", color = "black", size = 1)+
labs(title = "30 min: Induced", subtitle = "Dynamic")+
xlab("")+ 
ylab("log2(RPKM+1)") +
theme(legend.position="none")+
scale_x_continuous(breaks=ii) + 
theme(axis.text.x = element_text(angle = 45, hjust = 1,size=8),axis.text.y = element_text(size=8)) +
ylim(0,10)
print(plot1b)
dev.off()

# 1 hr
up1 <- GeneSets_LightSim_Deg[[2]][GeneSets_LightSim_Deg[[2]]$log2FoldChange < 0,]
down1 <- GeneSets_LightSim_Deg[[2]][GeneSets_LightSim_Deg[[2]]$log2FoldChange > 0,]


filtUp1 <- df[df$Rows %in% rownames(up1),]
filtUp1 <- droplevels(filtUp1)

filtDown1 <- df[df$Rows %in% rownames(down1),]
filtDown1 <- droplevels(filtDown1)


filtUp1 <- filtUp1[!(filtUp1$Rows %in% filtUp30$Rows),]
filtDown1 <- filtDown1[!(filtDown1$Rows %in% filtDown30$Rows),]


ii <- c(0,0.5,1,2,3,4,5,6)
pdf("visualizations/MultiLine_1hr_INDUCED.pdf",width=3,height=4,useDingbats=FALSE)
plot2a <- ggplot(filtUp1, aes(variable, log, color=Labels)) + 
geom_line(alpha=0.5) +
theme_classic()+
#scale_color_manual(values=rep("grey",length(filtUp30$Labels)))+
stat_summary(fun.y = "mean", geom = "line", color = "black", size = 1)+
labs(title = "1hr: Induced", subtitle = "Dynamic")+
xlab("")+ 
ylab("log2(RPKM+1)") +
theme(legend.position="none")+
scale_x_continuous(breaks=ii) + 
theme(axis.text.x = element_text(angle = 45, hjust = 1,size=8),axis.text.y = element_text(size=8)) +
ylim(0,10)
print(plot2a)
dev.off()


ii <- c(0,0.5,1,2,3,4,5,6)
pdf("visualizations/MultiLine_1hr_REDUCED.pdf",width=3,height=4,useDingbats=FALSE)
plot2b <- ggplot(filtDown1, aes(variable, log, color=Labels)) + 
geom_line(alpha=0.5) +
theme_classic()+
#scale_color_manual(values=rep("grey",length(filtUp30$Labels)))+
stat_summary(fun.y = "mean", geom = "line", color = "black", size = 1)+
labs(title = "1hr: Reduced", subtitle = "Dynamic")+
xlab("")+ 
ylab("log2(RPKM+1)") +
theme(legend.position="none")+
scale_x_continuous(breaks=ii) + 
theme(axis.text.x = element_text(angle = 45, hjust = 1,size=8),axis.text.y = element_text(size=8)) +
ylim(0,10)
print(plot2b)
dev.off()

# 3 hr
up3 <- GeneSets_LightSim_Deg[[3]][GeneSets_LightSim_Deg[[3]]$log2FoldChange < 0,]
down3 <- GeneSets_LightSim_Deg[[3]][GeneSets_LightSim_Deg[[3]]$log2FoldChange > 0,]


filtUp3 <- df[df$Rows %in% rownames(up3),]
filtUp3 <- droplevels(filtUp3)

filtDown3 <- df[df$Rows %in% rownames(down3),]
filtDown3 <- droplevels(filtDown3)

vec <- unique(c(as.character(filtUp30$Rows),as.character(filtUp1$Rows)))
filtUp3 <- filtUp3[!(filtUp3$Rows %in% vec),]
filtUp3 <- droplevels(filtUp3)

vec <- unique(c(as.character(filtDown30$Rows),as.character(filtDown1$Rows)))
filtDown3 <- filtDown3[!(filtDown3$Rows %in% vec),]
filtDown3 <- droplevels(filtDown3)


ii <- c(0,0.5,1,2,3,4,5,6)
pdf("visualizations/MultiLine_3hr_INDUCED.pdf",width=3,height=4,useDingbats=FALSE)
plot3a <- ggplot(filtUp3, aes(variable, log, color=Labels)) + 
geom_line(alpha=0.5) +
theme_classic()+
#scale_color_manual(values=rep("grey",length(filtUp30$Labels)))+
stat_summary(fun.y = "mean", geom = "line", color = "black", size = 1)+
labs(title = "3hr: Induced", subtitle = "Dynamic")+
xlab("")+ 
ylab("log2(RPKM+1)") +
theme(legend.position="none")+
scale_x_continuous(breaks=ii) + 
theme(axis.text.x = element_text(angle = 45, hjust = 1,size=8),axis.text.y = element_text(size=8)) +
ylim(0,10)
print(plot3a)
dev.off()


ii <- c(0,0.5,1,2,3,4,5,6)
pdf("visualizations/MultiLine_3hr_REDUCED.pdf",width=3,height=4,useDingbats=FALSE)
plot3b <- ggplot(filtDown1, aes(variable, log, color=Labels)) + 
geom_line(alpha=0.5) +
theme_classic()+
#scale_color_manual(values=rep("grey",length(filtUp30$Labels)))+
stat_summary(fun.y = "mean", geom = "line", color = "black", size = 1)+
labs(title = "3hr: Reduced", subtitle = "Dynamic")+
xlab("")+ 
ylab("log2(RPKM+1)") +
theme(legend.position="none")+
scale_x_continuous(breaks=ii) + 
theme(axis.text.x = element_text(angle = 45, hjust = 1,size=8),axis.text.y = element_text(size=8)) +
ylim(0,10)
print(plot3b)
dev.off()

# 6 hr
up6 <- GeneSets_LightSim_Deg[[4]][GeneSets_LightSim_Deg[[4]]$log2FoldChange < 0,]
down6 <- GeneSets_LightSim_Deg[[4]][GeneSets_LightSim_Deg[[4]]$log2FoldChange > 0,]

tmpUp6 <- GeneSets_LightSim_Deg[[5]][GeneSets_LightSim_Deg[[5]]$log2FoldChange < 0,]
tmpDown6 <- GeneSets_LightSim_Deg[[5]][GeneSets_LightSim_Deg[[5]]$log2FoldChange > 0,]


filtUp6 <- df[df$Rows %in% rownames(up6),]
filtUp6 <- droplevels(filtUp6)

filtDown6 <- df[df$Rows %in% rownames(down6),]
filtDown6 <- droplevels(filtDown6)

vec <- unique(c(as.character(filtUp30$Rows),as.character(filtUp1$Rows),as.character(filtUp3$Rows),as.character(rownames(tmpUp6))))
filtUp6 <- filtUp6[!(filtUp6$Rows %in% vec),]
filtUp6 <- droplevels(filtUp6)

vec <- unique(c(as.character(filtDown30$Rows),as.character(filtDown1$Rows),as.character(filtDown3$Rows),as.character(rownames(tmpDown6))))
filtDown6 <- filtDown6[!(filtDown6$Rows %in% vec),]
filtDown6 <- droplevels(filtDown6)

ii <- c(0,0.5,1,2,3,4,5,6)
pdf("visualizations/MultiLine_6hr_INDUCED.pdf",width=3,height=4,useDingbats=FALSE)
plot4a <- ggplot(filtUp6, aes(variable, log, color=Labels)) + 
geom_line(alpha=0.5) +
theme_classic()+
#scale_color_manual(values=rep("grey",length(filtUp30$Labels)))+
stat_summary(fun.y = "mean", geom = "line", color = "black", size = 1)+
labs(title = "6hr : Induced", subtitle = "Dynamic")+
xlab("")+ 
ylab("log2(RPKM+1)") +
theme(legend.position="none")+
scale_x_continuous(breaks=ii) + 
theme(axis.text.x = element_text(angle = 45, hjust = 1,size=8),axis.text.y = element_text(size=8)) +
ylim(0,10)
print(plot4a)
dev.off()


ii <- c(0,0.5,1,2,3,4,5,6)
pdf("visualizations/MultiLine_6hr_REDUCED.pdf",width=3,height=4,useDingbats=FALSE)
plot4b <- ggplot(filtDown6, aes(variable, log, color=Labels)) + 
geom_line(alpha=0.5) +
theme_classic()+
#scale_color_manual(values=rep("grey",length(filtUp30$Labels)))+
stat_summary(fun.y = "mean", geom = "line", color = "black", size = 1)+
labs(title = "6hr : Reduced", subtitle = "Dynamic")+
xlab("")+ 
ylab("log2(RPKM+1)") +
theme(legend.position="none")+
scale_x_continuous(breaks=ii) + 
theme(axis.text.x = element_text(angle = 45, hjust = 1,size=8),axis.text.y = element_text(size=8)) +
ylim(0,10)
print(plot4b)
dev.off()

plot <- plot_grid(plot1a,plot2a,plot3a,plot4a,labels=c("A","B","C","D"), ncol = 2,nrow = 2,align = "h")
save_plot("visualizations/MultiLine_COWPLOT_INDUCED.pdf", plot, ncol = 2,base_height=6,base_width=2)

plot <- plot_grid(plot1b,plot2b,plot3b,plot4b,labels=c("A","B","C","D"), ncol = 2,nrow = 2,align = "h")
save_plot("visualizations/MultiLine_COWPLOT_REDUCED.pdf", plot, ncol = 2,base_height=6,base_width=2)

