suppressPackageStartupMessages(library(org.Mm.eg.db))
suppressPackageStartupMessages(library(annotate))
suppressPackageStartupMessages(library(Rsubread))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(TxDb.Mmusculus.UCSC.mm10.knownGene))

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
genes <- transcriptsBy(txdb, "gene")
prom <- promoters(genes, upstream=1000,downstream=100)

genesDF <- as.data.frame(genes)
genesDF$Gene <- getSYMBOL(genesDF$group_name, data='org.Mm.eg')

result <- genesDF %>% 
             group_by(Gene) %>%
             filter(width == max(width)) %>%
             arrange(Gene) %>%
             as.data.frame()

promoters <- as.data.frame(prom)
promoters <- promoters[promoters$tx_id %in% result$tx_id,]

promoters$Gene <- getSYMBOL(promoters$group_name, data='org.Mm.eg')

promoters <- promoters %>%
             arrange(Gene) %>%
             as.data.frame() %>%
             select(-group,-group_name,-tx_id)

write.table(promoters,"Promoter_All_Genes_1kb.bed",sep="\t",row.names=F,col.names=F,quote=F)

finalgenes = read.table("inputs/Increase_halfhr_total.txt", header = T)
tmp <- promoters[promoters$Gene %in% rownames(finalgenes),] %>%
             select(-width,-tx_name) %>%
             arrange(seqnames)
write.table(tmp,"min30_1kb_Prom_Induced.bed",sep="\t",row.names=F,col.names=F,quote=F)

# 
finalgenes = read.table("inputs/Increase_1hr_unique.txt", header = T)
tmp <- promoters[promoters$Gene %in% finalgenes$Gene,]%>%
             select(-width,-tx_name) %>%
             arrange(seqnames)
write.table(tmp,"hr01_1kb_Prom_Induced.bed",sep="\t",row.names=F,col.names=F,quote=F)


# 
finalgenes = read.table("inputs/Increase_3hr_unique.txt", header = T)
tmp <- promoters[promoters$Gene %in% finalgenes$Gene,]%>%
             select(-width,-tx_name) %>%
             arrange(seqnames)
write.table(tmp,"hr03_1kb_Prom_Induced.bed",sep="\t",row.names=F,col.names=F,quote=F)

# 
finalgenes = read.table("inputs/Increase_6hr_unique.txt", header = T)
tmp <- promoters[promoters$Gene %in% finalgenes$Gene,]%>%
             select(-width,-tx_name) %>%
             arrange(seqnames)
write.table(tmp,"hr06_1kb_Prom_Induced.bed",sep="\t",row.names=F,col.names=F,quote=F)

# 
finalgenes = read.table("inputs/Decrease_1hr_unique.txt", header = T)
tmp <- promoters[promoters$Gene %in% finalgenes$Gene,]%>%
             select(-width,-tx_name) %>%
             arrange(seqnames)
write.table(tmp,"hr01_1kb_Prom_Reduced.bed",sep="\t",row.names=F,col.names=F,quote=F)

# 
finalgenes = read.table("inputs/Decrease_3hr_unique.txt", header = T)
tmp <- promoters[promoters$Gene %in% finalgenes$Gene,]%>%
             select(-width,-tx_name) %>%
             arrange(seqnames)
write.table(tmp,"hr03_1kb_Prom_Reduced.bed",sep="\t",row.names=F,col.names=F,quote=F)

# 
finalgenes = read.table("inputs/Decrease_6hr_unique.txt", header = T)
tmp <- promoters[promoters$Gene %in% finalgenes$Gene,]%>%
             select(-width,-tx_name) %>%
             arrange(seqnames)
write.table(tmp,"hr06_1kb_Prom_Reduced.bed",sep="\t",row.names=F,col.names=F,quote=F)