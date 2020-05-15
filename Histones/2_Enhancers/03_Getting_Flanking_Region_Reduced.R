suppressPackageStartupMessages(library(biomaRt))
suppressPackageStartupMessages(library(regioneR))
suppressPackageStartupMessages(library(BSgenome.Mmusculus.UCSC.mm10))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(igraph))

dir.create("Reduced_Flanking_Region/")

library(httr)
library(devtools)
set_config(
  use_proxy(url="proxy.swmed.edu", port=3128,username="S157784",password="Arianna2109")
)

# Getting the mart for mus musculus
mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
att <- listAttributes(mart)

# Getting the flanking region for 1 hr Ligth Stimulation Genes
genes = read.table("inputs/Dge_Time/Decrease_1hr_unique.txt", header = T)

transcripts <- getBM(attributes=c("chromosome_name", "start_position", "end_position", "ensembl_gene_id","strand", "gene_biotype", "mgi_symbol"),
      filters = "mgi_symbol", 
      values = genes$Gene, 
      uniqueRows = TRUE, 
        mart =mart)
transcripts$TSS = ifelse(transcripts$strand ==1, transcripts$start_position, transcripts$end_position)

flanking_50kb = GRanges(transcripts$mgi_symbol,
                           seqnames = Rle(paste0('chr', transcripts$chromosome_name)),
                           ranges = IRanges(start = transcripts$TSS - 50000, end = transcripts$TSS + 50000),
                           strand = Rle(rep("*", nrow(transcripts))))

flanking_100kb = GRanges(transcripts$mgi_symbol,
                           seqnames = Rle(paste0('chr', transcripts$chromosome_name)),
                           ranges = IRanges(start = transcripts$TSS - 100000, end = transcripts$TSS + 100000),
                           strand = Rle(rep("*", nrow(transcripts))))

#Change the chromosomes names from Ensembl (1,2,3...) to UCSC (chr1, chr2, ch3...) so they match the chromosome names in the BSgenome
seqlevelsStyle(flanking_50kb) <- "UCSC"
seqlevelsStyle(flanking_100kb) <- "UCSC"

df <- data.frame(seqnames=seqnames(flanking_50kb),
starts=start(flanking_50kb)-1,
ends=end(flanking_50kb),
gene=transcripts$mgi_symbol)
write.table(df,"Reduced_Flanking_Region/hr01_50kb_Flanking.bed",sep="\t",row.names=F,col.names=F,quote=F)

df <- data.frame(seqnames=seqnames(flanking_100kb),
starts=start(flanking_100kb)-1,
ends=end(flanking_100kb),
gene=transcripts$mgi_symbol)
write.table(df,"Reduced_Flanking_Region/hr01_100kb_Flanking.bed",sep="\t",row.names=F,col.names=F,quote=F)


# Getting the flanking region for 3 hr Ligth Stimulation Genes
genes = read.table("inputs/Dge_Time/Decrease_3hr_unique.txt", header = T)

transcripts <- getBM(attributes=c("chromosome_name", "start_position", "end_position", "ensembl_gene_id","strand", "gene_biotype", "mgi_symbol"),
      filters = "mgi_symbol", 
      values = genes$Gene, 
      uniqueRows = TRUE, 
        mart =mart)
transcripts$TSS = ifelse(transcripts$strand ==1, transcripts$start_position, transcripts$end_position)

flanking_50kb = GRanges(transcripts$mgi_symbol,
                           seqnames = Rle(paste0('chr', transcripts$chromosome_name)),
                           ranges = IRanges(start = transcripts$TSS - 50000, end = transcripts$TSS + 50000),
                           strand = Rle(rep("*", nrow(transcripts))))

flanking_100kb = GRanges(transcripts$mgi_symbol,
                           seqnames = Rle(paste0('chr', transcripts$chromosome_name)),
                           ranges = IRanges(start = transcripts$TSS - 100000, end = transcripts$TSS + 100000),
                           strand = Rle(rep("*", nrow(transcripts))))

#Change the chromosomes names from Ensembl (1,2,3...) to UCSC (chr1, chr2, ch3...) so they match the chromosome names in the BSgenome
seqlevelsStyle(flanking_50kb) <- "UCSC"
seqlevelsStyle(flanking_100kb) <- "UCSC"

df <- data.frame(seqnames=seqnames(flanking_50kb),
starts=start(flanking_50kb)-1,
ends=end(flanking_50kb),
gene=transcripts$mgi_symbol)
write.table(df,"Reduced_Flanking_Region/hr03_50kb_Flanking.bed",sep="\t",row.names=F,col.names=F,quote=F)

df <- data.frame(seqnames=seqnames(flanking_100kb),
starts=start(flanking_100kb)-1,
ends=end(flanking_100kb),
gene=transcripts$mgi_symbol)
write.table(df,"Reduced_Flanking_Region/hr03_100kb_Flanking.bed",sep="\t",row.names=F,col.names=F,quote=F)


# Getting the flanking region for 6 hr Ligth Stimulation Genes
genes = read.table("inputs/Dge_Time/Decrease_6hr_unique.txt", header = T)

transcripts <- getBM(attributes=c("chromosome_name", "start_position", "end_position", "ensembl_gene_id","strand", "gene_biotype", "mgi_symbol"),
      filters = "mgi_symbol", 
      values = genes$Gene, 
      uniqueRows = TRUE, 
        mart =mart)
transcripts$TSS = ifelse(transcripts$strand ==1, transcripts$start_position, transcripts$end_position)

flanking_50kb = GRanges(transcripts$mgi_symbol,
                           seqnames = Rle(paste0('chr', transcripts$chromosome_name)),
                           ranges = IRanges(start = transcripts$TSS - 50000, end = transcripts$TSS + 50000),
                           strand = Rle(rep("*", nrow(transcripts))))

flanking_100kb = GRanges(transcripts$mgi_symbol,
                           seqnames = Rle(paste0('chr', transcripts$chromosome_name)),
                           ranges = IRanges(start = transcripts$TSS - 100000, end = transcripts$TSS + 100000),
                           strand = Rle(rep("*", nrow(transcripts))))

#Change the chromosomes names from Ensembl (1,2,3...) to UCSC (chr1, chr2, ch3...) so they match the chromosome names in the BSgenome
seqlevelsStyle(flanking_50kb) <- "UCSC"
seqlevelsStyle(flanking_100kb) <- "UCSC"

df <- data.frame(seqnames=seqnames(flanking_50kb),
starts=start(flanking_50kb)-1,
ends=end(flanking_50kb),
gene=transcripts$mgi_symbol)
write.table(df,"Reduced_Flanking_Region/hr06_50kb_Flanking.bed",sep="\t",row.names=F,col.names=F,quote=F)

df <- data.frame(seqnames=seqnames(flanking_100kb),
starts=start(flanking_100kb)-1,
ends=end(flanking_100kb),
gene=transcripts$mgi_symbol)
write.table(df,"Reduced_Flanking_Region/hr06_100kb_Flanking.bed",sep="\t",row.names=F,col.names=F,quote=F)

