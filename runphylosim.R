#!/usr/bin/Rscript
library(phylosim)
library(multicore)
library(phangorn)
PSIM_FAST=TRUE
basetree <- rcoal(13,c("SE001","SE002","SE003","SE004","SE005","SE006","SE007","SE008","SE009","SE010","SE011","SE012","SE013"))
p <- WAG()
s <- AminoAcidSequence(length=200, processes=list(list(p)))
sampleStates(s)
d<- ContinuousDeletor(rate = 0.25, max.length = 10, dist = expression(rnorm(1, mean = 5, sd = 3)))
i<- ContinuousInsertor(rate = 0.25, max.length = 10, dist = expression(rnorm(1, mean = 5, sd = 3)))
i$templateSeq <- AminoAcidSequence(length = 10, processes = list(list(p)))
attachProcess(s,d)
attachProcess(s,i)
setDeletionTolerance(s, d, 0.08 + c(1/2:51, 1/51:2))
setInsertionTolerance(s, i, 0.08 + c(1/2:51, 1/51:2))


sim.replicate <- function(i) {
      name <- paste("replication_", i, sep = "")
      clearStates(s)
      sampleStates(s)
      sim <- Simulate(PhyloSim(name = name, root.seq = s, phylo = t), quiet = FALSE)
      if (i < 10) saveAlignment(sim, file = paste("gene00", i, ".fas", sep = ""),skip.internal=TRUE)
      else if (100 > i & i >= 10) saveAlignment(sim, file = paste("gene0", i, ".fas", sep = ""),skip.internal=TRUE)
      else if (1000 > i & i >= 100) saveAlignment(sim, file = paste("gene", i, ".fas", sep = ""),skip.internal=TRUE)
      return(sim)
      }



rep.size <- 6
lower <- 1
upper <- lower + rep.size - 1
t <- rNNI(basetree,1)
res.objects <- mclapply(lower:upper, sim.replicate)
print(res.objects)

lower <- upper + 1
upper <- lower + rep.size - 1
t <- rNNI(basetree,5)
res.objects <- mclapply(lower:upper, sim.replicate)
print(res.objects)

lower <- upper + 1
upper <- lower + rep.size - 1
t <- rNNI(basetree,5)
res.objects <- mclapply(lower:upper, sim.replicate)
print(res.objects)

lower <- upper + 1
upper <- lower + rep.size - 1
t <- rNNI(basetree,5)
res.objects <- mclapply(lower:upper, sim.replicate)
print(res.objects)

lower <- upper + 1
upper <- lower + rep.size - 1
t <- rNNI(basetree,5)
res.objects <- mclapply(lower:upper, sim.replicate)
print(res.objects)
