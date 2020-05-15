library(reshape2)
library(WGCNA)
library(tidyverse)

tab <- read.table("motifs_100kb.txt",header=T,sep="\t")
df <- reshape2::dcast(tab,Gene ~ Class,value.var="FDR",fun.aggregate = mean)
rownames(df) <- df$Gene
df$Gene <- NULL
df[is.na(df)] <- 1

order <- c("min30","hr01","hr03","hr06")
df <- df[,match(order,colnames(df))]

df <- -log10(df)


df1 <- df %>% rownames_to_column("TF") %>% arrange(desc(min30),desc(hr01),desc(hr03),desc(hr06)) %>% column_to_rownames("TF")

# Colors
low0 <- rgb(246, 248, 251, maxColorValue = 255)
low1 <- rgb(201, 230, 234, maxColorValue = 255)
low2 <- rgb(89, 158, 193, maxColorValue = 255)
mid <- rgb(67, 121, 180, maxColorValue = 255)
high1 <- rgb(65, 90, 158, maxColorValue = 255)
high2 <- rgb(52, 52, 106, maxColorValue = 255)
high3 <- rgb(15, 15, 25, maxColorValue = 255)


pdf("Motif_100kb_Cascade.pdf",width=4,height=10)
labeledHeatmap(Matrix = df1, xLabels = colnames(df1), yLabels = rownames(df1), 
	colorLabels =FALSE,colors=colorRampPalette(c(low0, low1, low2, mid, high1))(50),
	setStdMargins = FALSE, cex.text = 0.5)
dev.off()