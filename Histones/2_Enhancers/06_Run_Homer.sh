ls *_Induced_100kb_H3k27ac_ForEnrich.bed | parallel --progress --eta -j 20 'findMotifsGenome.pl {} mm10 {.}_HOMER -size 500 -mask -len 8,20,25 -p 15 -S 30 -cpg'

ls *_Reduced_100kb_H3k27ac_ForEnrich.bed | parallel --progress --eta -j 20 'findMotifsGenome.pl {} mm10 {.}_HOMER -size 500 -mask -len 8,20,25 -p 15 -S 30 -cpg'
