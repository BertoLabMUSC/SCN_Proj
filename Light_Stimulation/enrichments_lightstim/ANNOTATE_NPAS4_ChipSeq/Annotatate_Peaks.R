library("Cairo")
library("ChIPseeker")
library("clusterProfiler")
library("TxDb.Mmusculus.UCSC.mm10.knownGene")
library(org.Mm.eg.db)
library(annotate)
library(tidyverse)


txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
files <- list.files(".", pattern = ".bed", full.names = TRUE)
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
peakAnno <- lapply(files, annotatePeak, TxDb=txdb, tssRegion=c(-3000, 3000), verbose=FALSE)

Npas4_Cell <- as.data.frame(peakAnno[[1]])
Npas4_Tk <- as.data.frame(peakAnno[[2]])


Npas4_Cell$Gene <- getSYMBOL(Npas4_Cell$geneId, data='org.Mm.eg')
Npas4_Tk$Gene <- getSYMBOL(Npas4_Tk$geneId, data='org.Mm.eg')


Npas4_Cell <- Npas4_Cell %>%
				select(Gene,annotation)

Npas4_Tk <- Npas4_Tk %>%
				select(Gene,annotation)

openxlsx::write.xlsx(Npas4_Cell, file = "Npas4_Cell_Annotation.xlsx", colNames = TRUE, borders = "columns",sheetName="Cell")
openxlsx::write.xlsx(Npas4_Tk, file = "Npas4_Tk_Annotation.xlsx", colNames = TRUE, borders = "columns",sheetName="TK")