library(phangorn)
treesets <- function(taxa,n,boost,target){
    # outputs sets of three trees which are equally distant from each other
    # taxa is number of taxa
    # n is unweighted RF distance between trees
    # boost is a speedup (set = to n)
    # target is number of sets to create
    # branch lengths are drawn from Gamma(2.5,40)
    l <- list()
    t <- rtree(taxa)
    e <- length(t$edge.length)
    t$edge.length <- rep(1.0, times=e)
    while(length(l) <= target){
        trees <- rNNI(t,n/2+boost,10)
        for(i in 1:10){
            if(dist.topo(t,trees[[i]])==n){
                moretrees <- rNNI(trees[[i]],n/2+boost,20)
                for(j in 1:20){
                    if(dist.topo(t,moretrees[[j]])==n){
                        if(dist.topo(trees[[i]],moretrees[[j]])==n){
                            t$edge.length <- rgamma(e,shape=2.5,scale=40)
                            trees[[i]]$edge.length <- rgamma(e,shape=2.5,scale=40)
                            moretrees[[j]]$edge.length <- rgamma(e,shape=2.5,scale=40)
                            vec <- c(write.tree(t),write.tree(trees[[i]]),write.tree(moretrees[[j]]))
                            l[[length(l)+1]] <- vec
                        }
                    }
                }
            }
        }
    }
    return(l[1:target])
}
