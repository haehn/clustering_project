#!/usr/bin/Rscript

library(phangorn)
options (echo = FALSE)

args <- commandArgs(trailingOnly = TRUE)

tree <- read.tree(args[1])
method <- args[2]
moves <- args[3]

permute <- function (tree, method, moves) {
    if (method == "nni") {
        return(rNNI(tree,moves,1))
    }
    else if (nethod == "spr") {
        return(rSPR(tree,moves,1))
    }
}

out <- permute (tree, method, moves)

par(mfrow=c(1,2))
plot(tree)
plot(out)
