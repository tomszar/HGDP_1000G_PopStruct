# HGDP and 1000G population structure analysis    

This repository contains the results and the procedures from population genetic analyses done on the HGDP and 1000G reference samples.
Different analyses were done for both datasets together, therefore, we have included different pipelines and results. 
The merging of the HGDP and 1000G was done [elsewhere](https://tomszar.github.io/HGDP_1000G_Merge/).

## QC procedure

The QC procedure applied to the genotype data can be seen [here](https://nbviewer.jupyter.org/github/tomszar/HGDP_1000G_PopStruct/blob/master/Code/2018-05-CleanSNPs.ipynb).
In summary we:

1. Removed founders, that is, individuals with at least one parent in the dataset, and retained only autosomal chromosomes
2. Removed SNPs with missing call rates higher than 0.1
3. Removed SNPs with minor allele frequencies below 0.05
4. Removed SNPs with hardy-weinberg equilibrium p-values less than 1e-50
5. Removed samples with missing call rates higher than 0.1
6. Removed one arbitrary individual from any pairwise comparison with a pihat higher than 0.25 from an IBD estimation after LD prune. A pihat value of 0.25 is second degree relative.

After the QC procedure we split our dataset by chromosome for the phasing-IBD pipeline.
This dataset contains 3291 people, and 458912 variants.
Additionally, we generated a second dataset LD pruned, containing 3291 people, and 155533 variants, to be used for the PCA and Admixture analyses.

## Admixture and PCA

We ran an Admixture analysis and computed a PCA on the genotype data. 
We used an LD pruned dataset to ran both analyses.
The code for the Admixture analysis is [here](https://github.com/tomszar/HGDP_1000G_PopStruct/blob/master/Code/AdmixtureRun.sh), while the PCA is contained in the [QC procedure code](https://nbviewer.jupyter.org/github/tomszar/HGDP_1000G_PopStruct/blob/master/Code/2018-05-CleanSNPs.ipynb))

## Refined IBD

Another approach was to ran an IBD analysis on the genotype data. 
To do that, we first phased our data. 
The script can be seen [here](https://github.com/tomszar/HGDP_1000G_PopStruct/blob/master/Code/PhasingGenos.sh).
To keep consistency, we phased the HGDP and 1000G together, although the 1000 Genomes data can be found already phased.
