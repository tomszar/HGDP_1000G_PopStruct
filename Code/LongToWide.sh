#!/bin/bash
#Navigating to DataBases folder
cd ../DataBases/Genotypes
awk '
NR>1 \
{snps[$2]; clust[$3]; val[$2,$3]=$6} \
END {
	printf "%s", "rsid"; for (i in snps) printf "\t%s", i; print "";
for (j in clust) {printf "%s", j;
for (i in snps) printf "\t%s", val[i,j]; print ""} }' CleanGenos_strict.frq.strat > CleanGenos_strict_wide.frq.strat
