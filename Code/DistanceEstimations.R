####LIBRARIES####
library(tidyverse)
library(adegenet)
library(data.table) #for fread function
source("DataFrame2Matrix.R")

####DATABASES####
setwd('..')
path <- getwd()
setwd(paste(path, "/Results/RefinedIBD", sep = ""))

cm <- read_delim("cMmatrix.txt", delim = " ", col_names = FALSE, col_types = cols(
  X1 = col_character(),
  X2 = col_character(),
  X3 = col_double()
))

setwd(paste(path, "/Results/Counts", sep = ""))
allelecounts <- fread("hgdp1000ghg19_founders_geno01_maf_hwe_mind01_rel_pruned.frq.strat.wide", sep = ",",
                      header = TRUE, data.table = FALSE, index = "Pop")
rownames(allelecounts) <- allelecounts$Pop
allelecounts           <- allelecounts[,-c(1,2)]

setwd(paste(path, "/DataBases/PopInfo", sep = ""))
pops_hgdp <- read_delim("SampleInformation.txt", delim = "\t")[,c(1,3,5)]
colnames(pops_hgdp) <- c("ID", "pop", "super_pop")

pops_100g <- read_delim("integrated_call_samples_v3.20130502.ALL.panel", delim = "\t", col_types = cols(
  X1 = col_character(),
  X2 = col_character(),
  X3 = col_character()
))[,c(1:3)]
colnames(pops_100g)[1] <- c("ID")

#Changing the hgdp id names
for (i in 1:dim(pops_hgdp)[1] ) {
  nint   <- nchar(pops_hgdp[i,1])
  nzeros <- paste(integer(5-nint), collapse="")
  pops_hgdp[i,1] <- paste("HGDP", nzeros, pops_hgdp[i,1], sep = "")
}
pops <- bind_rows(pops_100g, pops_hgdp)

#Changing the HGDP super_pop
pops$super_pop[pops$super_pop == "AFRICA"]      <- "AFR"
pops$super_pop[pops$super_pop == "AMERICA"]     <- "AMR"
pops$super_pop[pops$super_pop == "EAST_ASIA"]   <- "EAS"
pops$super_pop[pops$super_pop == "EUROPE"]      <- "EUR"
pops$super_pop[pops$super_pop == "MIDDLE_EAST"] <- "MDE"
pops$super_pop[pops$super_pop == "OCEANIA"]     <- "OCE"
pops$super_pop[pops$super_pop == "CENTRAL_SOUTH_ASIA"] <- "SAS"

#Changing to factors
pops <- pops %>% mutate(pop=as.factor(pop), super_pop=as.factor(super_pop))

#counts to genpop
mygenpop <- genpop(allelecounts, ploidy=as.integer(2), type="codom")
rm(allelecounts)

####NEI DISTANCE####
#Getting distance
neidist <- dist.genpop(mygenpop)
rm(mygenpop)
#Saving matrix
setwd(paste(path, "/Results/Distances", sep = ""))
write.csv(as.matrix(neidist), file = "Nei_dist.csv", quote = FALSE)

####IBD DISTANCE####
myMatrix  <- data.frame2matrix(cm, 'X1', 'X2', 'X3')
rm(cm)
ind <- lower.tri(myMatrix)
myMatrix[ind]  <- t(myMatrix)[ind] 
diag(myMatrix) <- 0
rm(ind)

#Getting an empty matrix
npop     <- length(attr(neidist, "Labels"))
popnames <- attr(neidist, "Labels")
pop_ibd_matrix <- matrix(0, ncol = npop, nrow = npop, dimnames = list(popnames, popnames ))

for (i in 1:npop ){
  popname1 <- popnames[i]
  d <- pops %>% filter(pop == popname1) %>% select(ID)
  d <- d$ID
  for (j in 1:npop ) {
    popname2 <- popnames[j]
    e <- pops %>% filter(pop == popname2) %>% select(ID)
    e <- e$ID
    t <- sum(myMatrix[ d[d %in% rownames(myMatrix)],  e[e %in% rownames(myMatrix)]]) / ( length(d) * length(e) )
    pop_ibd_matrix[popname1, popname2] <- t
  }
}

rm(myMatrix)
diag(pop_ibd_matrix) <- 0
#Getting sqrt values
pop_ibd_matrix <- sqrt(pop_ibd_matrix)
#Getting similarity matrix
pop_ibd_matrix <- max(pop_ibd_matrix) - pop_ibd_matrix

#Saving matrix
diag(pop_ibd_matrix) <- 0
setwd(paste(path, "/Results/Distances", sep = ""))
write.csv(as.matrix(pop_ibd_matrix), file = "IBD_dist.csv", quote = FALSE)