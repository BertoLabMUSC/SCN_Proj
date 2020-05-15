# Analysis with Enhancers H3k27ac
mkdir "Reduced_Enhancer"

# hr01
bedtools intersect -a Reduced_Flanking_Region/hr01_100kb_Flanking.bed -b H3k27ac_Enhancers/H3k27ac_Consensus_Enhancers.bed -c > Reduced_Enhancer/hr01_Reduced_100kb_H3k27ac.bed

# hr03
bedtools intersect -a Reduced_Flanking_Region/hr03_100kb_Flanking.bed -b H3k27ac_Enhancers/H3k27ac_Consensus_Enhancers.bed -c > Reduced_Enhancer/hr03_Reduced_100kb_H3k27ac.bed

# hr06
bedtools intersect -a Reduced_Flanking_Region/hr06_100kb_Flanking.bed -b H3k27ac_Enhancers/H3k27ac_Consensus_Enhancers.bed -c > Reduced_Enhancer/hr06_Reduced_100kb_H3k27ac.bed

# For enrichment

# hr01
bedtools intersect -a H3k27ac_Enhancers/H3k27ac_Consensus_Enhancers.bed -b Reduced_Flanking_Region/hr01_100kb_Flanking.bed -c | awk -v OFS="\t" '{if($4 > 0 ){print}}' > Reduced_Enhancer/hr01_Reduced_100kb_H3k27ac_ForEnrich.bed

# hr03
bedtools intersect -a H3k27ac_Enhancers/H3k27ac_Consensus_Enhancers.bed -b Reduced_Flanking_Region/hr03_100kb_Flanking.bed -c | awk -v OFS="\t" '{if($4 > 0 ){print}}' > Reduced_Enhancer/hr03_Reduced_100kb_H3k27ac_ForEnrich.bed

# hr06
bedtools intersect -a H3k27ac_Enhancers/H3k27ac_Consensus_Enhancers.bed -b Reduced_Flanking_Region/hr06_100kb_Flanking.bed -c | awk -v OFS="\t" '{if($4 > 0 ){print}}' > Reduced_Enhancer/hr06_Reduced_100kb_H3k27ac_ForEnrich.bed
