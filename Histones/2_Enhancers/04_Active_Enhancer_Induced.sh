# Analysis with Enhancers H3k27ac
mkdir "Induced_Enhancer"

# 30min
bedtools intersect -a Induced_Flanking_Region/min30_100kb_Flanking.bed -b H3k27ac_Enhancers/H3k27ac_Consensus_Enhancers.bed -c > Induced_Enhancer/min30_Induced_100kb_H3k27ac.bed

# hr01
bedtools intersect -a Induced_Flanking_Region/hr01_100kb_Flanking.bed -b H3k27ac_Enhancers/H3k27ac_Consensus_Enhancers.bed -c > Induced_Enhancer/hr01_Induced_100kb_H3k27ac.bed

# hr03
bedtools intersect -a Induced_Flanking_Region/hr03_100kb_Flanking.bed -b H3k27ac_Enhancers/H3k27ac_Consensus_Enhancers.bed -c > Induced_Enhancer/hr03_Induced_100kb_H3k27ac.bed

# hr06
bedtools intersect -a Induced_Flanking_Region/hr06_100kb_Flanking.bed -b H3k27ac_Enhancers/H3k27ac_Consensus_Enhancers.bed -c > Induced_Enhancer/hr06_Induced_100kb_H3k27ac.bed

# For enrichment
bedtools intersect -a H3k27ac_Enhancers/H3k27ac_Consensus_Enhancers.bed -b Induced_Flanking_Region/min30_100kb_Flanking.bed -c | awk -v OFS="\t" '{if($4 > 0 ){print}}' > Induced_Enhancer/min30_Induced_100kb_H3k27ac_ForEnrich.bed

# hr01
bedtools intersect -a H3k27ac_Enhancers/H3k27ac_Consensus_Enhancers.bed -b Induced_Flanking_Region/hr01_100kb_Flanking.bed -c | awk -v OFS="\t" '{if($4 > 0 ){print}}' > Induced_Enhancer/hr01_Induced_100kb_H3k27ac_ForEnrich.bed

# hr03
bedtools intersect -a H3k27ac_Enhancers/H3k27ac_Consensus_Enhancers.bed -b Induced_Flanking_Region/hr03_100kb_Flanking.bed -c | awk -v OFS="\t" '{if($4 > 0 ){print}}' > Induced_Enhancer/hr03_Induced_100kb_H3k27ac_ForEnrich.bed

# hr06
bedtools intersect -a H3k27ac_Enhancers/H3k27ac_Consensus_Enhancers.bed -b Induced_Flanking_Region/hr06_100kb_Flanking.bed -c | awk -v OFS="\t" '{if($4 > 0 ){print}}' > Induced_Enhancer/hr06_Induced_100kb_H3k27ac_ForEnrich.bed
