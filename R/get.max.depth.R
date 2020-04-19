get.max.depth <-
function(phy){
    ht=heights.phylo(phy)
	N=Ntip(phy)
	Tmax=ht$start[N+1]
    Tmax
}
