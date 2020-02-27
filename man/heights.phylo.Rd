\name{heights.phylo}
\alias{heights.phylo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
heights.phylo(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x) 
{
    phy = x
    phy <- reorder(phy, "postorder")
    n <- length(phy$tip.label)
    n.node <- phy$Nnode
    xx <- numeric(n + n.node)
    for (i in nrow(phy$edge):1) xx[phy$edge[i, 2]] <- xx[phy$edge[i, 
        1]] + phy$edge.length[i]
    root = ifelse(is.null(phy$root.edge), 0, phy$root.edge)
    labs = c(phy$tip.label, phy$node.label)
    depth = max(xx)
    tt = depth - xx
    idx = 1:length(tt)
    dd = phy$edge.length[idx]
    mm = match(1:length(tt), c(phy$edge[, 2], Ntip(phy) + 1))
    dd = c(phy$edge.length, root)[mm]
    ss = tt + dd
    res = cbind(ss, tt)
    rownames(res) = idx
    colnames(res) = c("start", "end")
    res = data.frame(res)
    res
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
