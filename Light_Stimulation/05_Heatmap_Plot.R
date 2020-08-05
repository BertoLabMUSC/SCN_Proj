library(ggplot2)
library(ggpubr)
library(ggforce)
library(reshape2)
library(tidyverse)
library(pheatmap)
library(ComplexHeatmap)

dir.create("HEATMAPS/")

tmpData <- read.table("inputs/RPKM_LighSim_EXON.txt")
ord=c("CT17_1","CT17_1s","CT17_2","CT17_2s","Min30_1","Min30_1s","Min30_2","Min30_2s","Hour1_1","Hour1_2","Hour3_1","Hour3_2","Hour6_1","Hour6_2","CT23_1","CT23_2")
tmpData <- tmpData[,match(ord,colnames(tmpData))]

# Heatmaps for Increased Genes
pdf("visualizations/MegaHeatmap_Increased.pdf",width=3,height=5,useDingbats=FALSE)
up1 <- read.table("inputs/Increase_halfhr_total.txt",header=T)
up2 <- read.table("inputs/Increase_1hr_unique.txt",header=T)
up3 <- read.table("inputs/Increase_3hr_unique.txt",header=T)
up4 <- read.table("inputs/Increase_6hr_unique.txt",header=T)

df <- rbind(up1,up2,up3,up4)
filtUp <- log2(tmpData[rownames(tmpData) %in% df$Gene,]+1)


mat_scaled = t(apply(filtUp, 1, scale))
colnames(mat_scaled) <- colnames(filtUp)

Heatmap(mat_scaled, 
            name = "z-score", 
            #row_km = 4,
            cluster_columns = FALSE,
            border = TRUE,
            show_row_names = FALSE)

dev.off()



# Heatmaps for Decreased Genes
pdf("visualizations/MegaHeatmap_Decreased.pdf",width=3,height=5,useDingbats=FALSE)
#up1 <- read.table("inputs/Decrease_halfhr_total.txt",header=T)
up2 <- read.table("inputs/Decrease_1hr_unique.txt",header=T)
up3 <- read.table("inputs/Decrease_3hr_unique.txt",header=T)
up4 <- read.table("inputs/Decrease_6hr_unique.txt",header=T)

df <- rbind(up2,up3,up4)
filtUp <- log2(tmpData[rownames(tmpData) %in% df$Gene,]+1)
filtUp$Min30 <- NULL

mat_scaled = t(apply(filtUp, 1, scale))
colnames(mat_scaled) <- colnames(filtUp)

mat_scaled <- mat_scaled %>% 
                  as.data.frame() %>% 
                  rownames_to_column("TF") %>% 
                  arrange(CT17,Hour1,Hour3,Hour6) %>% 
                  column_to_rownames("TF")



Heatmap(as.matrix(mat_scaled), 
            name = "z-score", 
            #row_km = 4,
            cluster_columns = FALSE,
            border = TRUE,
            show_row_names = FALSE)

dev.off()

# Heatmaps for Decreased Genes with 30 min
pdf("visualizations/MegaHeatmap_Decreased_WithMin30.pdf",width=3,height=5,useDingbats=FALSE)
up1 <- read.table("inputs/Decrease_halfhr_total.txt",header=T)
up2 <- read.table("inputs/Decrease_1hr_unique.txt",header=T)
up3 <- read.table("inputs/Decrease_3hr_unique.txt",header=T)
up4 <- read.table("inputs/Decrease_6hr_unique.txt",header=T)

df <- rbind(up1,up2,up3,up4)
filtUp <- log2(tmpData[rownames(tmpData) %in% df$Gene,]+1)
#filtUp$Min30 <- NULL

mat_scaled = t(apply(filtUp, 1, scale))
colnames(mat_scaled) <- colnames(filtUp)

mat_scaled <- mat_scaled %>% 
                  as.data.frame() %>% 
                  rownames_to_column("TF") %>% 
                  arrange(desc(CT17),desc(Hour1),desc(Hour3),desc(Hour6)) %>% 
                  column_to_rownames("TF")

Heatmap(as.matrix(mat_scaled), 
            name = "z-score", 
            #row_km = 4,
            cluster_columns = FALSE,
            border = TRUE,
            show_row_names = FALSE)

dev.off()



