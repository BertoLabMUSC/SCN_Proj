# CLUSTERING DATA | ALL CELLS
This script performs data clustering on single nuclei RNA-seq UMI counts dataset in multiple steps such as filtering, normalization, scaling using Seurat pipeline (Seurat v3.0). Towards the end, gene markers per cluster are also determined using seurat functions. *sessionInfo* is provided towards the end of the script.



#### Data Clustering

###### Load Counts Data and Create Seurat Object

```R
######################### Load collapsed and combibed data
load("SCN_DATA_COLLAPSED_COMBINED.RData")
dim(allData)
## total 144,806 cells

######################### Create Seurat object without any filtering
seuObj <- CreateSeuratObject(counts = allData, project = "SCN")

######################### Create meta data
metaData <- as.data.frame(seuObj@meta.data)
colnames(metaData) <- c("Library", "nUMI", "nGenes")

######################### Create conditions
metaData$Condition <- rep("NA", nrow(metaData))
metaData[grepl("^D", metaData[,"Library"]),"Condition"] <- "Dark"
metaData[grepl("^L", metaData[,"Library"]),"Condition"] <- "Light"

######################### Create batches
metaData$Batch <- rep("0", nrow(metaData))
metaData[grepl("D0|L0", metaData[,"Library"]),"Batch"] <- "1"
metaData[grepl("D1|D2|D3|D4|L1|L2|L3|L4", metaData[,"Library"]),"Batch"] <- "2"

######################### Add library
all.lib <- metaData$Library
names(all.lib) <- row.names(metaData)

######################### Add condition
all.cond <- metaData$Condition
names(all.cond) <- row.names(metaData)

######################### Add batch
all.batch <- metaData$Batch
names(all.batch) <- row.names(metaData)

######################### Update seurat object
seuObj <- AddMetaData(object = seuObj, metadata = all.lib, col.name = "Library")
seuObj <- AddMetaData(object = seuObj, metadata = all.cond, col.name = "Condition")
seuObj <- AddMetaData(object = seuObj, metadata = all.batch, col.name = "Batch")

######################### Percent mitochondrial content
mito.genes <- grep(pattern = "^mt-", x = rownames(x = seuObj@assays$RNA@data), value = TRUE)
percent.mito <- Matrix::colSums(seuObj@assays$RNA@data[mito.genes, ])/Matrix::colSums(seuObj@assays$RNA@data)
seuObj <- AddMetaData(object = seuObj, metadata = percent.mito, col.name = "pMito")

######################### Save RData
save(seuObj, file = "SCN_SEURAT_DATA.RData")
```



###### Data QC

```R
######################### Violin plots colored by Libraries
plot_nU <- VlnPlot(object = seuObj, features = "nCount_RNA", pt.size = 0, group.by = "Library")
plot_nG <- VlnPlot(object = seuObj, features = "nFeature_RNA", pt.size = 0, group.by = "Library")
plot_pM <- VlnPlot(object = seuObj, features = "pMito", pt.size = 0, group.by = "Library")
plotQC1 <- grid.arrange(plot_nU, plot_nG, plot_pM, ncol = 3)
ggsave(paste(seuObj@project.name, "_QC_1_LIBRARY.pdf", sep = ""), plot = plotQC1, width = 16, height = 4, units = "in", dpi = 300)

######################### Violin plots colored by Condition
plot_nU <- VlnPlot(object = seuObj, features = "nCount_RNA", pt.size = 0, group.by = "Condition")
plot_nG <- VlnPlot(object = seuObj, features = "nFeature_RNA", pt.size = 0, group.by = "Condition")
plot_pM <- VlnPlot(object = seuObj, features = "pMito", pt.size = 0, group.by = "Condition")
plotQC1 <- grid.arrange(plot_nU, plot_nG, plot_pM, ncol = 3)
ggsave(paste(seuObj@project.name, "_QC_1_CONDITION.pdf", sep = ""), plot = plotQC1, width = 12, height = 4, units = "in", dpi = 300)

######################### Violin plots colored by Batch
plot_nU <- VlnPlot(object = seuObj, features = "nCount_RNA", pt.size = 0, group.by = "Batch")
plot_nG <- VlnPlot(object = seuObj, features = "nFeature_RNA", pt.size = 0, group.by = "Batch")
plot_pM <- VlnPlot(object = seuObj, features = "pMito", pt.size = 0, group.by = "Batch")
plotQC1 <- grid.arrange(plot_nU, plot_nG, plot_pM, ncol = 3)
ggsave(paste(seuObj@project.name, "_QC_1_BATCH.pdf", sep = ""), plot = plotQC1, width = 12, height = 4, units = "in", dpi = 300)

######################### Scatter plots colored by Libraries
plot_nUpM <- FeatureScatter(object = seuObj, feature1 = "nCount_RNA", feature2 = "pMito", cex.use = 1, group.by = "Library")
plot_nUnG <- FeatureScatter(object = seuObj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", cex.use = 1, group.by = "Library")
plotQC2 <- grid.arrange(plot_nUpM, plot_nUnG, ncol = 2)
ggsave(paste(seuObj@project.name, "_QC_2_LIBRARY.pdf", sep = ""), plot = plotQC2, width = 10, height = 4, units = "in", dpi = 300)

######################### Scatter plots colored by Condition
plot_nUpM <- FeatureScatter(object = seuObj, feature1 = "nCount_RNA", feature2 = "pMito", cex.use = 1, group.by = "Condition")
plot_nUnG <- FeatureScatter(object = seuObj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", cex.use = 1, group.by = "Condition")
plotQC2 <- grid.arrange(plot_nUpM, plot_nUnG, ncol = 2)
ggsave(paste(seuObj@project.name, "_QC_2_CONDITION.pdf", sep = ""), plot = plotQC2, width = 10, height = 4, units = "in", dpi = 300)

######################### Scatter plots colored by Batch
plot_nUpM <- FeatureScatter(object = seuObj, feature1 = "nCount_RNA", feature2 = "pMito", cex.use = 1, group.by = "Batch")
plot_nUnG <- FeatureScatter(object = seuObj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", cex.use = 1, group.by = "Batch")
plotQC2 <- grid.arrange(plot_nUpM, plot_nUnG, ncol = 2)
ggsave(paste(seuObj@project.name, "_QC_2_BATCH.pdf", sep = ""), plot = plotQC2, width = 10, height = 4, units = "in", dpi = 300)

```



###### Filtering Data

```R
######################### Data filtering cutoffs
nUlo <- -Inf ## lower cutoff for number of UMI
nUhi <- 10000 ## higher cutoff for number of UMI
nGlo <- -Inf ## lower cutoff for number of genes
nGhi <- 4000 ## higher cutoff for number of genes
pMlo <- -Inf ## lower cutoff for percent mito content
pMhi <- 0.1 #(10%) ## higher cutoff for percent mito content

######################### Scatter plots colored by Libraries (post-filtering)
plotseuObjAll_UM <- ggplot(seuObj@meta.data, aes(x = nCount_RNA, y = pMito, colour = Library)) + geom_point(shape = 20, size = 2, alpha = 1) + labs(title = paste("ALL_DATA", round(cor(seuObj@meta.data$nCount_RNA, seuObj@meta.data$pMito), 2), sep = ", "), x = "Number of UMIs", y = "Percent Mito") + annotate("rect", xmin = nUlo, xmax = nUhi, ymin = -Inf, ymax = Inf, alpha = 0.1, fill = "blue") + annotate("rect", xmin = -Inf, xmax = Inf, ymin = pMlo, ymax = pMhi, alpha = 0.1, fill = "green") + theme_bw()
plotseuObjAll_UG <- ggplot(seuObj@meta.data, aes(x = nCount_RNA, y = nFeature_RNA, colour = Library)) + geom_point(shape = 20, size = 2, alpha = 1) + labs(title = paste("ALL_DATA", round(cor(seuObj@meta.data$nCount_RNA, seuObj@meta.data$nFeature_RNA), 2), sep = ", "), x = "Number of UMIs", y = "Number of Genes") + annotate("rect", xmin = nUlo, xmax = nUhi, ymin = -Inf, ymax = Inf, alpha = 0.1, fill = "red") + annotate("rect", xmin = -Inf, xmax = Inf, ymin = nGlo, ymax = nGhi, alpha = 0.1, fill = "green") + theme_bw()
plotQCseuObjAll2 <- grid.arrange(plotseuObjAll_UM, plotseuObjAll_UG, ncol = 2)
ggsave(paste(seuObj@project.name, "_FILTER.pdf", sep = ""), plot = plotQCseuObjAll2, width = 14, height = 6, units = "in", dpi = 150)

######################### Applying filters
seuObjFiltTemp <- SubsetData(object = seuObj, subset.name = "nCount_RNA", low.threshold = nUlo, high.threshold = nUhi)
seuObjFiltTemp2 <- SubsetData(object = seuObjFiltTemp, subset.name = "nFeature_RNA", low.threshold = nGlo, high.threshold = nGhi)
seuObjFilt <- SubsetData(object = seuObjFiltTemp2, subset.name = "pMito", low.threshold = pMlo, high.threshold = pMhi)
dim(seuObjFilt@assays$RNA@data)
## 144,582 filtered cells out of total 144,806 cells

######################### Genes from chr X Y amd M
xgenes <- scan("GENCODE_vM17_MM10_ChrX_GENES.txt", what = "", sep = "\n")
ygenes <- scan("GENCODE_vM17_MM10_ChrY_GENES.txt", what = "", sep = "\n")
mgenes <- scan("GENCODE_vM17_MM10_ChrM_GENES.txt", what = "", sep = "\n")
genes2remove <- sort(c(mgenes, xgenes, ygenes))

######################### Create the list of genes and cells to retain after filtering
keepGenes <- unique(sort(setdiff(row.names(seuObjFilt@assays$RNA@data), genes2remove)))
keepCells <- colnames(seuObjFilt@assays$RNA@data)

######################### Filter UMI data using genes and cells to retain
scnDataFilt <- seuObjFilt@assays$RNA@counts[keepGenes, keepCells]

######################### Create seurat object with the filtered raw (non-normalized data)
seuObjFilt <- CreateSeuratObject(counts = scnDataFilt, project = "SCN_FILT")

######################### Update meta data for pMito
pMitoFilt <- seuObj@meta.data[keepCells, "pMito"]
names(pMitoFilt) <- row.names(seuObj@meta.data[keepCells,])
seuObjFilt <- AddMetaData(object = seuObj, metadata = pMitoFilt, col.name = "pMito")

######################### Update meta data for library
libFilt <- seuObj@meta.data[keepCells, "Library"]
names(libFilt) <- row.names(seuObj@meta.data[keepCells,])
seuObjFilt <- AddMetaData(object = seuObj, metadata = libFilt, col.name = "Library")

######################### Update meta data for genotype
condFilt <- seuObj@meta.data[keepCells, "Condition"]
names(condFilt) <- row.names(seuObj@meta.data[keepCells,])
seuObjFilt <- AddMetaData(object = seuObj, metadata = condFilt, col.name = "Condition")

table(seuObjFilt@meta.data$Library)
#    D0    D1    D2    D3    D4    L0    L1    L2    L3    L4
# 14131 11768 15972 18569 13397  9056 11154 15099 21405 14255

table(seuObjFilt@meta.data$Condition)
#  Dark Light
# 73731 70851

######################### Save RData
save(seuObjFilt, keepCells, keepGenes, file = "SCN_SEURAT_FILT.RData")
```



###### Normaliza and Scale Data

```R
######################### LogNormalize the data using scaling factor
seuObjFiltNorm <- NormalizeData(object = seuObjFilt, normalization.method = "LogNormalize", scale.factor = 10000)

######################### Identify ~2000 variable genes
seuObjFiltNorm <- FindVariableFeatures(object = seuObjFiltNorm)

######################### Scale the data
seuObjFiltNorm <- ScaleData(object = seuObjFiltNorm, vars.to.regress = c("nCount_RNA", "pMito", "Batch"), model.use = "linear")

######################### Save RData files
save(seuObjFiltNorm, file = "SCN_SEURAT_FILT_NORM.RData")
```



###### PCA Analysis

```R
######################### PCA using identified variable genes on regressed (scaled) data
seuObjFiltNormPC <- RunPCA(object = seuObjFiltNorm, pc.genes = seuObjFiltNorm@assays$RNA@var.features, do.print = TRUE, pcs.print = 1:5, genes.print = 5, pcs.compute = 50)

######################### Examine and visualize PCA results a few different ways
pdf("SCN_SEURAT_FILT_NORM_PCA_1.pdf", width = 9, height = 21)
VizDimLoadings(object = seuObjFiltNormPC, pcs.use = 1:30)
dev.off()

pdf("SCN_SEURAT_FILT_NORM_2.pdf", width = 9, height = 6)
DimPlot(object = seuObjFiltNormPC, dim.1 = 1, dim.2 = 2)
dev.off()

seuObjFiltNormPC <- ProjectDim(object = seuObjFiltNormPC, reduction = "pca")

pdf("SCN_SEURAT_FILT_NORM_PCA_3.pdf", width = 9, height = 21)
DimHeatmap(object = seuObjFiltNormPC, dims = 1:12, cells = 500, balanced = TRUE)
dev.off()

pdf("SCN_SEURAT_FILT_NORM_PCA_4.pdf", width = 6, height = 4)
ElbowPlot(object = seuObjFiltNormPC, ndims = 50, reduction = "pca")
dev.off()

######################### Determine statistically significant principal components
seuObjFiltNormPC <- JackStraw(object = seuObjFiltNormPC, reduction = "pca", dims = 50, num.replicate = 100)
seuObjFiltNormPC <- ScoreJackStraw(object = seuObjFiltNormPC, dims = 1:50)

pdf("SCN_SEURAT_FILT_NORM_PCA_5.pdf", width = 12, height = 8)
JackStrawPlot(object = seuObjFiltNormPC, dims = 1:50, reduction = "pca")
dev.off()

######################### Save RData
save(seuObjFiltNormPC, file = "SCN_SEURAT_FILT_NORM_PCA.RData")
```



###### Data Clustering

```R
######################### Find neighbors & clusters
numPCs <- 50
seuObjFiltNormPCClust <- FindNeighbors(object = seuObjFiltNormPC, reduction = "pca", dims = 1:numPCs)
seuObjFiltNormPCClust <- FindClusters(object = seuObjFiltNormPCClust, save.SNN = TRUE)

######################### Run UMAP
seuObjFiltNormPCClust <- RunUMAP(object = seuObjFiltNormPCClust, reduction = "pca", dims = 1:numPCs)

pdf("SCN_SEURAT_FILT_NORM_PCA_CLUST_UMAP.pdf", width = 7, height = 6)
DimPlot(object = seuObjFiltNormPCClust, reduction = "umap", label = TRUE, label.size = 4, repel = T, pt.size = 1)
dev.off()

######################### Save RData
save(seuObjFiltNormPCClust, file = "SCN_SEURAT_FILT_NORM_PCA_CLUST.RData")
```



#### Cluster Markers

```R
######################### Find cluster markers
seuObjFiltNormPCClustDEG <- FindAllMarkers(object = seuObjFiltNormPCClust, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)

######################### Save cluster markers and RData
write.table(seuObjFiltNormPCClustDEG, "SCN_SEURAT_FILT_NORM_PCA_CLUST_DEG_TABLE.txt", row.names = T, col.names = T, quote = F, sep = "\t")
save(seuObjFiltNormPCClustDEG, file = "SCN_SEURAT_FILT_NORM_PCA_CLUST_DEG.RData")
```



#### R Session Info
```R
sessionInfo()

R version 3.5.1 (2018-07-02)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Red Hat Enterprise Linux Server 7.4 (Maipo)

Matrix products: default
BLAS/LAPACK: /cm/shared/apps/intel/compilers_and_libraries/2017.6.256/linux/mkl/lib/intel64_lin/libmkl_gf_lp64.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] RColorBrewer_1.1-2 reticulate_1.12    ggrepel_0.8.0      reshape2_1.4.3    
[5] gridExtra_2.3      ggplot2_3.1.1      Matrix_1.2-14      dplyr_0.8.0.1     
[9] Seurat_3.0.0.9000 

loaded via a namespace (and not attached):
 [1] httr_1.4.0          tidyr_0.8.3         jsonlite_1.6       
 [4] viridisLite_0.3.0   splines_3.5.1       lsei_1.2-0         
 [7] R.utils_2.8.0       gtools_3.8.1        Rdpack_0.11-0      
[10] assertthat_0.2.1    globals_0.12.4      pillar_1.3.1       
[13] lattice_0.20-35     glue_1.3.1          digest_0.6.18      
[16] SDMTools_1.1-221.1  colorspace_1.4-1    cowplot_0.9.4      
[19] htmltools_0.3.6     R.oo_1.22.0         plyr_1.8.4         
[22] pkgconfig_2.0.2     bibtex_0.4.2        tsne_0.1-3         
[25] listenv_0.7.0       purrr_0.3.2         scales_1.0.0       
[28] RANN_2.6.1          gdata_2.18.0        Rtsne_0.15         
[31] tibble_2.1.1        withr_2.1.2         ROCR_1.0-7         
[34] pbapply_1.4-0       lazyeval_0.2.2      survival_2.42-3    
[37] magrittr_1.5        crayon_1.3.4        R.methodsS3_1.7.1  
[40] future_1.12.0       nlme_3.1-137        MASS_7.3-50        
[43] gplots_3.0.1.1      ica_1.0-2           tools_3.5.1        
[46] fitdistrplus_1.0-14 data.table_1.12.2   gbRd_0.4-11        
[49] stringr_1.4.0       plotly_4.9.0        munsell_0.5.0      
[52] cluster_2.0.7-1     irlba_2.3.3         compiler_3.5.1     
[55] rsvd_1.0.0          caTools_1.17.1.2    rlang_0.3.1        
[58] grid_3.5.1          ggridges_0.5.1      htmlwidgets_1.3    
[61] igraph_1.2.4.1      bitops_1.0-6        npsurv_0.4-0       
[64] gtable_0.3.0        codetools_0.2-15    R6_2.4.0           
[67] zoo_1.8-5           future.apply_1.2.0  KernSmooth_2.23-15 
[70] metap_1.1           ape_5.3             stringi_1.4.3      
[73] parallel_3.5.1      Rcpp_1.0.0          png_0.1-7          
[76] tidyselect_0.2.5    lmtest_0.9-37

```

