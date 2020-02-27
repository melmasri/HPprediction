\name{rEta.copheneticFast}
\alias{rEta.copheneticFast}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
rEta.copheneticFast(eta.old, tree, tree.ht, pdist.old, no0, i, sZ, Z, ywU, eta_sd = 0.01, a, b, nr, ind, sparse, dummyMat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eta.old}{
%%     ~~Describe \code{eta.old} here~~
}
  \item{tree}{
%%     ~~Describe \code{tree} here~~
}
  \item{tree.ht}{
%%     ~~Describe \code{tree.ht} here~~
}
  \item{pdist.old}{
%%     ~~Describe \code{pdist.old} here~~
}
  \item{no0}{
%%     ~~Describe \code{no0} here~~
}
  \item{i}{
%%     ~~Describe \code{i} here~~
}
  \item{sZ}{
%%     ~~Describe \code{sZ} here~~
}
  \item{Z}{
%%     ~~Describe \code{Z} here~~
}
  \item{ywU}{
%%     ~~Describe \code{ywU} here~~
}
  \item{eta_sd}{
%%     ~~Describe \code{eta_sd} here~~
}
  \item{a}{
%%     ~~Describe \code{a} here~~
}
  \item{b}{
%%     ~~Describe \code{b} here~~
}
  \item{nr}{
%%     ~~Describe \code{nr} here~~
}
  \item{ind}{
%%     ~~Describe \code{ind} here~~
}
  \item{sparse}{
%%     ~~Describe \code{sparse} here~~
}
  \item{dummyMat}{
%%     ~~Describe \code{dummyMat} here~~
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
function (eta.old, tree, tree.ht, pdist.old, no0, i, sZ, Z, ywU, 
    eta_sd = 0.01, a, b, nr, ind, sparse, dummyMat) 
{
    change = FALSE
    eta.prop = eta.old + eta_sd * rnorm(1)
    dist = cophFast(eb.phylo(tree, tree.ht, eta.prop), a, b, 
        nr, dummyMat)
    pdist.new = dist[i, ] \%*\% sZ
    if (sparse) 
        pdist.new = pdist.new@x
    if (length(no0)) {
        if (length(no0) == sum(Z[i, ])) 
            likeli = -Inf
        else likeli = sum(((log(pdist.new) - log(pdist.old)) * 
            Z[i, ])[-no0]) - sum((ywU * (pdist.new - pdist.old))[-no0])
    }
    else {
        likeli = sum((log(pdist.new) - log(pdist.old)) * Z[i, 
            ]) - sum(ywU * (pdist.new - pdist.old))
    }
    if (!is.na(likeli) && runif(1) <= min(1, exp(likeli))) {
        eta.old = eta.prop
        pdist.old = c(pdist.new)
        change = TRUE
    }
    list(eta = eta.old, dist = pdist.old, change = change)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
