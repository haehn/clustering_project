#!/usr/bin/Rscript
library(phylosim)
library(multicore)
library(phangorn)
PSIM_FAST=TRUE

make.basetree <- function(ntax) {
      vect = vector()
      for (i in 1:ntax) {
            if (i < 10) vect[i] <- paste("SE00",i,sep="")
            else if (100 > i & i >= 10) vect[i] <- paste("SE0",i,sep="")
            else if (1000 > i & i >= 100) vect[i] <- paste("SE",i,sep="")
      }
      return (rcoal(ntax,vect))
}

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

sim.classes <- function(n.reps,rep.size,n.nnis,basetree) {
      lower <- 1
      for (i in 1:n.reps) {
              upper <- lower + rep.size - 1
              print (lower)
              print(upper)
              t <- rNNI(basetree,n.nnis)
              #res.objects <- mclapply(lower:upper, sim.replicate)
              #print(res.objects)
              lower <- upper + 1

      }
}

basetree <- make.basetree(5)
#basetree <- rcoal(8,c("SE001","SE002","SE003","SE004","SE005","SE006","SE007","SE008"))
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



#sim.classes(2,4,4,basetree)

rep.size <- 10
lower <- 1
upper <- lower + rep.size - 1
t <- rNNI(basetree,1)
res.objects <- mclapply(lower:upper, sim.replicate,mc.preschedule=TRUE, mc.cores=8)
print(res.objects)

# lower <- upper + 1
# upper <- lower + rep.size - 1
# t <- rNNI(basetree,5)
# res.objects <- mclapply(lower:upper, sim.replicate)
# print(res.objects)

# lower <- upper + 1
# upper <- lower + rep.size - 1
# t <- rNNI(basetree,5)
# res.objects <- mclapply(lower:upper, sim.replicate)
# print(res.objects)

# lower <- upper + 1
# upper <- lower + rep.size - 1
# t <- rNNI(basetree,5)
# res.objects <- mclapply(lower:upper, sim.replicate)
# print(res.objects)

# # lower <- upper + 1
# # upper <- lower + rep.size - 1
# # t <- rNNI(basetree,5)
# # res.objects <- mclapply(lower:upper, sim.replicate)
# # print(res.objects)
