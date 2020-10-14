get.max.depth <-function(phy){
    if(is.phylo(phy)){
        ht=heights.phylo(phy)
        N=Ntip(phy)
        Tmax=ht$start[N+1]
        Tmax
    }else{
        max(phy)
    }
}
