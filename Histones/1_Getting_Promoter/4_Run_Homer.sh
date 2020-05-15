ls *_Induced.bed | parallel --progress --eta -j 4 'findMotifsGenome.pl {} mm10 {.}_HOMER -mask -len 8,20,25 -p 15 -S 30 -cpg'

ls *_Reduced.bed | parallel --progress --eta -j 4 'findMotifsGenome.pl {} mm10 {.}_HOMER -mask -len 8,20,25 -p 15 -S 30 -cpg'

