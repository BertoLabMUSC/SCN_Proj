# COMBINE DATA
This SCN single nuclei RNA-seq data was generated for two conditions (i) Dark and (ii) Light. Each condition has 5 biological replicates. Single nuclei RNA-seq libraries were generated for each sample using 10X Genomics v2 microfluidics platform. And each library was sequenced more than once (technical replicates) to increase the reads depth. Purpose of this script is to collapse the technical replicates for each library for each condition by selecting common cells across technical replicates followed by combining all libraries from all conditions into a single UMI counts matrix for downstream processing. The code chunks may need to be re-run changing the library/condition name in order to process all libraries for all conditions.



#### COLLAPSE TECHINCAL REPLICATES PER LIBRARY/SAMPLE

###### Read Reference Genes ID/Names List (Gencode vM17)
Raw UMI counts matrix for each sequencing run of every library is not identical due to many reasons. The goal of this part is to match all technical replicates using a reference list of gene names/symbols from Gencode vM17 reference annotation.
```R
refGenes <- read.table("gencode.vM17.protein_coding_gene_id_symbol.txt", sep = "\t", header = T)
refGenesSymbol <- as.data.frame(unique(sort(refGenes$GeneSymbol)))
colnames(refGenesSymbol) <- "genes"
row.names(refGenesSymbol) <- refGenesSymbol$genes
print(dim(refGenesSymbol))
```


###### Identify, extract and collapse common cells for each library
For each library, common cell barcodes are first identified and a plot showing common and unique cell barcodes across technical replicates is generated using *UpSetR* R package. Using the list of common cell barcodes identified across technical replicates of each library, data corresponding to common cell barcodes is extracted and collapsed using *aggregate* function. And, aggregated data frame for each library is stored as RData object.
```R
collapseCommonCells <- function(dataTR1, dataTR2, nameLIB)
	{
	print(paste("Processing:", nameLIB, sep = " "))
	print(paste("Data Technical Rep 1:", ncol(dataTR1), "Cells", sep = " "))
	print(paste("Data Technical Rep 2:", ncol(dataTR2), "Cells", sep = " "))

	######################### IDENTIFY COMMON CELLS
	libCells <- list(TR1 = colnames(dataTR1), TR2 = colnames(dataTR2))

	pdf(paste("SCN_10X", nameLIB, "Common_CellBarcodes.pdf", sep = "_"), width = 6, height = 4)
	libPlot <- upset(fromList(libCells), order.by = "freq")
	dev.off()

	libCombos <- Reduce(c, lapply(2:length(libCells), function(x) combn(1:length(libCells), x, simplify=FALSE) ))
	libIntersect <- lapply(libCombos, function(x) Reduce(intersect, libCells[x]) )

	commonCellsLib <- libIntersect[[max(length(libIntersect))]]
	print(paste("Common Cells:", length(commonCellsLib), sep = " "))

	write.table(commonCellsLib, paste("SCN_10X", nameLIB, "Common_CellBarcodes.txt", sep = "_"), row.names = F, col.names = F, quote = F, sep = "\t")


	######################### EXTRACT COMMON CELLS
	dataTR1common <- dataTR1[,colnames(dataTR1) %in% commonCellsLib]
	colnames(dataTR1common) <- paste(nameLIB, "TR1", colnames(dataTR1common), sep = "_")
	dataTR1common$genes <- row.names(dataTR1common)

	dataTR2common <- dataTR2[,colnames(dataTR2) %in% commonCellsLib]
	colnames(dataTR2common) <- paste(nameLIB, "TR2", colnames(dataTR2common), sep = "_")
	dataTR2common$genes <- row.names(dataTR2common)


	######################### COLLAPSE COMMON CELLS
	combinedData <- merge(dataTR1common, dataTR2common, by = "genes")
	row.names(combinedData) <- combinedData$genes
	combinedData$genes <- NULL

	combinedMeta <- as.data.frame(matrix(unlist(strsplit(colnames(combinedData), "_")), ncol = 3, byrow = TRUE))
	row.names(combinedMeta) <- colnames(combinedData)
	colnames(combinedMeta) <- c("Library", "Technical", "CellBarcode")

	newDataTemp <- as.data.frame(t(combinedData))
	newData <- merge(newDataTemp, combinedMeta, by = "row.names")
	row.names(newData) <- newData$Row.names
	newData$Row.names <- NULL
	cellBC <- list(newData$CellBarcode)
	newData$Library <- NULL
	newData$CellBarcode <- NULL
	newData$Technical <- NULL

	newDataAggrTemp <- aggregate(newData, by = cellBC, FUN = "sum")
	row.names(newDataAggrTemp) <- newDataAggrTemp$Group.1
	newDataAggrTemp$Group.1 <- NULL

	newDataAggr <- as.data.frame(t(newDataAggrTemp))
	colnames(newDataAggr) <- paste(nameLIB, colnames(newDataAggr), sep = "_")
	newDataAggr$genes <- row.names(newDataAggr)

	######################### SAVE COLLAPSED DATA
	save(newDataAggr, file = paste(nameLIB, "Collapsed_Data.RData", sep = "_"))
	print(paste("Collapse Complete:", nameLIB, sep = ""))
	}


for (lib in c("Lib01", "Lib02", "Lib03", "Lib04", "Lib05", "Lib06", "Lib07", "Lib08", "Lib09", "Lib10"))
	{
	print(paste("=========> Starting", lib, sep = ""))
	data1 <- read.table(<path_to_data_folder_1_lib>, header = T, sep = "\t", row.names = 1)
	print(dim(data1))
	data2 <- read.table(<path_to_data_folder_2_lib>, header = T, sep = "\t", row.names = 1)
	print(dim(data2))

	collapseCommonCells(data1, data2, lib)
	print(paste("=========> Done", lib, sep = ""))
	}

```


###### Read and update collapsed UMI count tables
Read collapsed UMI counts table for each library into a list and merge into a massive data frame with reference gene symbols. For the genes missing in the counts table, replace NAs with zeros.
```R
######################### Reference Genes
refGenes <- read.table("gencode.vM17.protein_coding_gene_id_symbol.txt", sep = "\t", header = T)
refGenesSymbol <- as.data.frame(unique(sort(refGenes$GeneSymbol)))
colnames(refGenesSymbol) <- "genes"
row.names(refGenesSymbol) <- refGenesSymbol$genes
print(dim(refGenesSymbol))

######################### Collapsed Data
print("Lib01")
liblist1 <- load("Lib01_Collapsed_Data.RData")
lib01.data <- newDataAggr
rm(list = liblist1)
rm(liblist1)
print(dim(lib01.data))
## Repeat for each library before next step

######################### List of collapsed Data with reference gene list
allDataTemp <- list(GENES = refGenesSymbol,
                    D0 = lib01.data,
                    D1 = lib02.data,
                    D2 = lib03.data,
                    D3 = lib04.data,
                    D4 = lib05.data,
                    L0 = lib06.data,
                    L1 = lib07.data,
                    L2 = lib08.data,
                    L3 = lib09.data,
                    L4 = lib10.data)

######################### Merge into single data frame
allData <- Reduce(function(x, y) { merge(x, y, all = TRUE, by = "genes") } , allDataTemp)
row.names(allData) <- allData$genes
allData$genes <- NULL
print(dim(allData))

######################### Replace NAs with Zeros
allData[is.na(allData)] <- 0
print(dim(allData))

######################### Save Collapsed and Combined Data
save(allData, file = "SCN_DATA_COLLAPSED_COMBINED.RData")
```



###### R Session Info
```
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
[1] gridExtra_2.3 ggplot2_3.1.1 UpSetR_1.3.3 

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.1       withr_2.1.2      assertthat_0.2.1 crayon_1.3.4    
 [5] dplyr_0.8.1      grid_3.5.1       plyr_1.8.4       R6_2.4.0        
 [9] gtable_0.3.0     magrittr_1.5     scales_1.0.0     pillar_1.4.0    
[13] rlang_0.3.4      lazyeval_0.2.2   glue_1.3.1       purrr_0.2.5     
[17] munsell_0.5.0    compiler_3.5.1   pkgconfig_2.0.2  colorspace_1.4-1
[21] tidyselect_0.2.5 tibble_2.1.1
```

