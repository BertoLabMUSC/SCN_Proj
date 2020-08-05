suppressPackageStartupMessages(library("Cairo"))
suppressPackageStartupMessages(library("ChIPseeker"))
suppressPackageStartupMessages(library("clusterProfiler"))
suppressPackageStartupMessages(library("TxDb.Mmusculus.UCSC.mm10.knownGene"))

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
files <- list.files(".", pattern = ".bed", full.names = FALSE)
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrix <- lapply(files, getTagMatrix, windows=promoter)
tagMatrixList <- setNames(tagMatrix,c("Consensus","Dark","Light"))
peakAnno <- lapply(files, annotatePeak, TxDb=txdb, tssRegion=c(-3000, 3000), verbose=FALSE)
peakAnnoList <- setNames(peakAnno,c("Consensus","Dark","Light"))

pdf("Pie_Chart_Consensus.pdf",width=7,height=7)
plotAnnoPie(peakAnno[[1]])
dev.off()

pdf("Pie_Chart_Dark.pdf",width=7,height=7)
plotAnnoPie(peakAnno[[2]])
dev.off()

pdf("Pie_Chart_Light.pdf",width=7,height=7)
plotAnnoPie(peakAnno[[3]])
dev.off()

pdf("BarChart_Annotation.pdf",width=6,height=4)
plotAnnoBar(peakAnnoList)
dev.off()



pdf("Distribution.pdf",width=5,height=5)
plotAvgProf(tagMatrixList, xlim=c(-3000, 3000), conf=0.95,resample=100, facet="row")
dev.off()

pdf("Profiles.pdf",width=5,height=7)
tagHeatmap(tagMatrixList, xlim=c(-3000, 3000), color=NULL)
dev.off()

