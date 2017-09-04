# Download mammal supertree to working directory
# Mammal Supertree
# Fritz et al. 2009 Ecology Letters
# 10.1111/j.1461-0248.2009.01307.x

require(ape)

if(!file.exists("example/mammal_supertree.tre")){
	require(fulltext)
	tree <- read.nexus(ft_get_si("10.1111/j.1461-0248.2009.01307.x", 1, "wiley"))[[1]]
	write.tree(tree, file="example/mammal_supertree.tre")
} else tree <- read.tree('example/mammal_supertree.tre')
