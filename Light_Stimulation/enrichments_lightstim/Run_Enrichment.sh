# Enrichment.r
# -g = Two columns table of input genes with specific association from your study
# -l = list of two columns tables with gene - disease association. E.g. Gene1 SYN
# -p = make a bubble chart with OR and -log10(FDR)
# -b = background (protein coding = 19776)
# -o = output label for statistics and viz
# -W/-H = width/height of the plot. 

mkdir enrichments_stats/

# Allen scRNA
Rscript Enrichment.r -g LightStim_PolyA_DEG.txt -l Npas4_Tk_GeneSets.RData -p -b 22000 -o enrichments_stats/Npas4_Tk -W 5 -H 3
Rscript Enrichment.r -g LightStim_PolyA_DEG.txt -l Npas4_Cell_GeneSets.RData -p -b 22000 -o enrichments_stats/Npas4_Cell_GeneSets -W 5 -H 3
Rscript Enrichment.r -g LightStim_PolyA_DEG.txt -l Greenberg_GeneSets.RData -p -b 22000 -o enrichments_stats/Greenberg_Targets -W 5 -H 3
Rscript Enrichment.r -g LightStim_PolyA_DEG.txt -l GeneSets_Enhancer.RData -p -b 22000 -o enrichments_stats/Enhancer -W 5 -H 3

