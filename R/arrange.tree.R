arrange.tree <-
function(phy){
    ## taken from .b.phylo from
    ## https://github.com/mwpennell/geiger-v2/blob/master/R/utilities-phylo.R#L1353-L1382
    ht=heights.phylo(phy)
	N=Ntip(phy)
	Tmax=ht$start[N+1]
	mm=match(1:nrow(ht), phy$edge[,2])
	ht$t1=Tmax-ht$end[phy$edge[mm,1]]
	ht$t2=ht$start-ht$end+ht$t1
    ht
}
