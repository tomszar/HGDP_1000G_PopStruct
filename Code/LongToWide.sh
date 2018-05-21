#!/bin/bash
#Navigating to DataBases folder
cd ../DataBases/Genotypes
awk '
NR>1 \
{snps[$2]; clust[$3]; a1[$2]=$4; a2[$2]=$5; val[$2,$3]=$7; total[$3]=$8} \
END {
	printf "%s,%s", "Pop", "Genomes"; for (i in snps) printf ",%s.%s,%s.%s", i,a1[i], i,a2[i]; print "";
for (j in clust) {printf "%s,%s", j, total[j];
for (i in snps) printf ",%s,%s", val[i,j], total[j] - val[i,j]; print ""} }' CleanGenos_strict.frq.strat > CleanGenos_strict_wide.frq.strat
