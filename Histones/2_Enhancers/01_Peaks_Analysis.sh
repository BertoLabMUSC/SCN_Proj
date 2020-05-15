# Consensus peakset
intersectBed -a inputs/Dark-H3k27ac_Broad_peaks.broadPeak -b inputs/Light-H3k27ac_Broad_peaks.broadPeak -f 0.2 -r | cut -f 1-3 > inputs/H3k27ac_Consensus.bed

# Intersect NPAS4
mkdir H3k27ac_Npas4_Intersections/
intersectBed -a inputs/H3k27ac_Consensus.bed -b inputs/Npas4_Tk_Narrow_peaks.bed -f 0.2 -r > H3k27ac_Npas4_Intersections/H3k27ac_Npas4_Consensus_Enhancers.bed

# Remove promoters
mkdir H3k27ac_Enhancers/
bedtools subtract -a inputs/H3k27ac_Consensus.bed -b inputs/Promoter_All_Genes_1kb.bed -f 0.2 -r > H3k27ac_Enhancers/H3k27ac_Consensus_Enhancers.bed

# Intersect promoters
mkdir H3k27ac_Promoters/
intersectBed -a inputs/H3k27ac_Consensus.bed -b inputs/Promoter_All_Genes_1kb.bed -f 0.2 -r > H3k27ac_Promoters/H3k27ac_Consensus_Promoter.bed
