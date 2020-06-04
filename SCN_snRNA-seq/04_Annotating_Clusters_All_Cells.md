# ANNOTATING CLUSTERS
This script covers a Fisher's exact test approach used to annotate and classify cells into different cell-types. Annotation is based on three published mouse single cell datasets.



#### Using Mouse Visual Cortex Dataset as Reference

Mouse visual cortex single cell RNA-seq dataset (Hrvatin et. al., 2017) was used to identify cluster-specific gene markers using Seurat (v3.0) pipeline retaining the original cell-barcode and cell-type association. Marker genes enriched for each cell-type were then used to perform Fisher's exact test to identify significant overlap of marker genes across our identified SCN clusters.

```R
######################### Read cluster markers
degtable <- list.files(path = "./../", pattern = 'SCN_SEURAT_FILT_NORM_PCA_CLUST_DEG_TABLE.txt')
tab <- read.table(paste("./../", degtable, sep = ""), sep = "\t", header = T)
tab <- tab[tab$p_val_adj <= 0.05 & tab$pct.1 >= 0.75,]
tab <- tab[c(7,8)]
tab$cluster <- as.factor(paste("Cluster_", sprintf("%02d", tab$cluster), sep=""))
colnames(tab) <- c("Gene","DEFINITION")
tab$DEFINITION <- as.factor(tab$DEFINITION)
Genes <- as.data.frame(table(tab$DEFINITION))

######################### Read reference cell-type markers
load("Hrvatin_VisualCortex_Data_QC_Filt_Norm_PCA_Clust_DEG.RData")

hrvatinAllSig <- hrvatin.markers[hrvatin.markers$p_val_adj <= 0.05,]
hrvatinAllSig2 <- hrvatinAllSig[hrvatinAllSig$pct.1 >= 0.5,]
hrvatinAllSig3 <- hrvatinAllSig2[c(7,6)]

vcGenes <- list(Astro = hrvatinAllSig3[hrvatinAllSig3$cluster == "Astro",],
                Endo_1 = hrvatinAllSig3[hrvatinAllSig3$cluster == "Endo_1",],
                Endo_2 = hrvatinAllSig3[hrvatinAllSig3$cluster == "Endo_2",],
                Macrophage = hrvatinAllSig3[hrvatinAllSig3$cluster == "Macrophage",],
                Micro_1 = hrvatinAllSig3[hrvatinAllSig3$cluster == "Micro_1",],
                Micro_2 = hrvatinAllSig3[hrvatinAllSig3$cluster == "Micro_2",],
                Neu_ExcL23 = hrvatinAllSig3[hrvatinAllSig3$cluster == "ExcL23",],
                Neu_ExcL4 = hrvatinAllSig3[hrvatinAllSig3$cluster == "ExcL4",],
                Neu_ExcL5_1 = hrvatinAllSig3[hrvatinAllSig3$cluster == "ExcL5_1",],
                Neu_ExcL5_2 = hrvatinAllSig3[hrvatinAllSig3$cluster == "ExcL5_2",],
                Neu_ExcL5_3 = hrvatinAllSig3[hrvatinAllSig3$cluster == "ExcL5_3",],
                Neu_ExcL6 = hrvatinAllSig3[hrvatinAllSig3$cluster == "ExcL6",],
                Neu_Int_Cck = hrvatinAllSig3[hrvatinAllSig3$cluster == "Int_Cck",],
                Neu_Int_Npy = hrvatinAllSig3[hrvatinAllSig3$cluster == "Int_Npy",],
                Neu_Int_Pv = hrvatinAllSig3[hrvatinAllSig3$cluster == "Int_Pv",],
                Neu_Int_Sst_1 = hrvatinAllSig3[hrvatinAllSig3$cluster == "Int_Sst_1",],
                Neu_Int_Sst_2 = hrvatinAllSig3[hrvatinAllSig3$cluster == "Int_Sst_2",],
                Neu_Int_Vip = hrvatinAllSig3[hrvatinAllSig3$cluster == "Int_Vip",],
                Olig_1 = hrvatinAllSig3[hrvatinAllSig3$cluster == "Olig_1",],
                Olig_2 = hrvatinAllSig3[hrvatinAllSig3$cluster == "Olig_2",],
                Olig_3 = hrvatinAllSig3[hrvatinAllSig3$cluster == "Olig_3",],
                Olig_4 = hrvatinAllSig3[hrvatinAllSig3$cluster == "Olig_4",],
                Olig_5 = hrvatinAllSig3[hrvatinAllSig3$cluster == "Olig_5",],
                Olig_6 = hrvatinAllSig3[hrvatinAllSig3$cluster == "Olig_6",],
                Olig_7 = hrvatinAllSig3[hrvatinAllSig3$cluster == "Olig_7",],
                OPC_1 = hrvatinAllSig3[hrvatinAllSig3$cluster == "OPC_1",],
                OPC_2 = hrvatinAllSig3[hrvatinAllSig3$cluster == "OPC_2",],
                Pericyte = hrvatinAllSig3[hrvatinAllSig3$cluster == "Pericyte",],
                SM_1 = hrvatinAllSig3[hrvatinAllSig3$cluster == "SM_1",],
                SM_2 = hrvatinAllSig3[hrvatinAllSig3$cluster == "SM_2",])

GeneSets <- vcGenes

######################### Re-arrange data
for(i in 1:length(GeneSets))
	{
	colnames(GeneSets[[i]])[1] <- "Gene"
	}

ln <- length(GeneSets)
cl <- length(Genes$Var1)
TEMP <- list()
INT <- list()
for (i in 1:ln)
	{
  TEMP[[i]] <- tab[tab$Gene %in% GeneSets[[i]]$Gene,]
  INT[[i]] <- as.data.frame(table(TEMP[[i]]$DEFINITION))
  }
names(INT) <- names(GeneSets)
names(TEMP) <- names(GeneSets)

######################### Create the matrix for each GeneSet
NROWS <- sapply(GeneSets, nrow)

#               GeneSet
#              +       -
#          +-------+-------+
#        + |   a   |   b   |
#  Module  +-------+-------+
#        - |   c   |   d   |
#          +-------+-------+

for (i in 1:length(INT))
  {
  INT[[i]]$b <- NROWS[[i]] - INT[[i]]$Freq
  INT[[i]]$c <- Genes$Freq - INT[[i]]$Freq
  INT[[i]]$d <- 9274 - Genes$Freq - nrow(GeneSets[[i]])
  }

######################### Function for Fisher's exact test
RunFisher <- function(row, alt = 'greater', cnf = 0.85) 
	{
  f <- fisher.test(matrix(row, nrow = 2), alternative = alt, conf.level = cnf)
  return(c(row,
           P_val = f$p.value,
           LogP = -log10(f$p.value),
           OR = f$estimate[[1]],
           OR_Low = f$conf.int[1],
           OR_Up = f$conf.int[2]))
	}

######################### Run Fisher's exact test
FisherMat=list()
for (i in 1:length(INT))
  {
  FisherMat[[i]] <- t(apply(INT[[i]][,2:5], 1, RunFisher))
  rownames(FisherMat[[i]]) <- INT[[i]]$Var1
  FisherMat[[i]] <- FisherMat[[i]][rownames(FisherMat[[i]]) != "grey",]
  }
names(FisherMat) <- names(INT)
save(FisherMat, TEMP, file = "SCN_SEURAT_FILT_NORM_PCA_CLUST_FISHER_OUTPUT_HRVATIN.RData")


######################### Arrange a matrix for P-val
tmp <- list()
FisherP <- matrix()
rowNames <- rownames(FisherMat[[1]])
colNames <- names(FisherMat)

for (i in 1:length(INT))
  {
  tmp[[i]] <- cbind(as.numeric(FisherMat[[i]][,5]))
  FisherP <- do.call(cbind,tmp)
  }
rownames(FisherP) <- rowNames
colnames(FisherP) <- colNames

######################### Arrange a matrix for OR
tmp <- list()
FisherOR <- matrix()
rowNames <- rownames(FisherMat[[1]])
colNames <- names(FisherMat)
for (i in 1:length(INT))
  {
  tmp[[i]] <- cbind(as.numeric(FisherMat[[i]][,6]))
  FisherOR <- do.call(cbind,tmp)
  }
rownames(FisherOR) <- rowNames
colnames(FisherOR) <- colNames

######################### Compute adjusted P-val
library(magrittr)
FisherAdj <- FisherP %>% as.matrix %>% as.vector %>% p.adjust(method='fdr') %>% matrix(ncol=ncol(FisherP))
rownames(FisherAdj) <- rowNames
colnames(FisherAdj) <- colNames

FisherAdj[FisherAdj > 0.05] <- 1
FisherOR[FisherOR < 1] <- 0

df <- -log10(FisherAdj)


######################### Plot
pdf("SCN_SEURAT_FILT_NORM_PCA_CLUST_FISHER_OUTPUT_HRVATIN_1.pdf",width=16, height=16, pointsize=15)
par(mar = c(11, 7, 2, 2))
LabelMat <- paste(signif(FisherOR, 2), "\n(",signif(FisherAdj, 1), ")", sep = "")
LabelMat[LabelMat == "0\n(1)"] <- NA
labeledHeatmap(Matrix = df,
               xLabels = colnames(df),
               yLabels = rownames(df),
               colorLabels =FALSE,
               colors=colorRampPalette(c("white", "red"))(50),
               textMatrix=LabelMat,
               setStdMargins = FALSE,
               cex.text = 0.5,
               xLabelsAngle = 90)
dev.off()


######################### Plot with different color scale
low0 <- rgb(246, 248, 251, maxColorValue = 255)
low1 <- rgb(201, 230, 234, maxColorValue = 255)
low2 <- rgb(89, 158, 193, maxColorValue = 255)
mid <- rgb(67, 121, 180, maxColorValue = 255)
high1 <- rgb(65, 90, 158, maxColorValue = 255)
high2 <- rgb(52, 52, 106, maxColorValue = 255)
high3 <- rgb(15, 15, 25, maxColorValue = 255)

pdf("SCN_SEURAT_FILT_NORM_PCA_CLUST_FISHER_OUTPUT_HRVATIN_2.pdf",width=12, height=16, pointsize=15)
par(mar = c(11, 7, 2, 2))
labeledHeatmap(Matrix = df,
              xLabels = colnames(df),
              yLabels = rownames(df),
              colorLabels =FALSE,
              colors=colorRampPalette(c(low0, low1, low2, mid, high1, high2, high3))(100),
              setStdMargins = FALSE,
              cex.text = 0.5,
              xLabelsAngle = 90)
dev.off()
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
[1] magrittr_1.5          dplyr_0.8.1           RColorBrewer_1.1-2   
[4] reshape2_1.4.3        ggplot2_3.1.1         cluster_2.0.7-1      
[7] WGCNA_1.68            fastcluster_1.1.25    dynamicTreeCut_1.63-1

loaded via a namespace (and not attached):
 [1] Biobase_2.42.0        bit64_0.9-7           splines_3.5.1        
 [4] foreach_1.4.4         Formula_1.2-3         assertthat_0.2.1     
 [7] stats4_3.5.1          latticeExtra_0.6-28   blob_1.1.1           
[10] fit.models_0.5-14     robustbase_0.93-3     impute_1.56.0        
[13] pillar_1.4.0          RSQLite_2.1.1         backports_1.1.4      
[16] lattice_0.20-35       glue_1.3.1            digest_0.6.19        
[19] checkmate_1.9.3       colorspace_1.4-1      htmltools_0.3.6      
[22] preprocessCore_1.44.0 Matrix_1.2-14         plyr_1.8.4           
[25] pcaPP_1.9-73          pkgconfig_2.0.2       purrr_0.2.5          
[28] GO.db_3.7.0           mvtnorm_1.0-8         scales_1.0.0         
[31] htmlTable_1.12        tibble_2.1.1          IRanges_2.16.0       
[34] withr_2.1.2           nnet_7.3-12           BiocGenerics_0.28.0  
[37] lazyeval_0.2.2        survival_2.42-6       crayon_1.3.4         
[40] memoise_1.1.0         doParallel_1.0.14     MASS_7.3-51          
[43] foreign_0.8-71        tools_3.5.1           data.table_1.12.2    
[46] matrixStats_0.54.0    stringr_1.4.0         S4Vectors_0.20.1     
[49] munsell_0.5.0         AnnotationDbi_1.44.0  compiler_3.5.1       
[52] rlang_0.3.4           grid_3.5.1            iterators_1.0.10     
[55] rstudioapi_0.10       htmlwidgets_1.3       robust_0.4-18        
[58] base64enc_0.1-3       gtable_0.3.0          codetools_0.2-15     
[61] DBI_1.0.0             rrcov_1.4-4           R6_2.4.0             
[64] gridExtra_2.3         knitr_1.20            bit_1.1-14           
[67] Hmisc_4.1-1           stringi_1.4.3         parallel_3.5.1       
[70] Rcpp_1.0.1            rpart_4.1-13          acepack_1.4.1        
[73] DEoptimR_1.0-8        tidyselect_0.2.5 

```

