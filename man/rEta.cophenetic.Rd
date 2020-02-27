\name{rEta.cophenetic}
\alias{rEta.cophenetic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
rEta.cophenetic(eta.old, tree, tree.ht, pdist.old, pd0, i, Z, ywU, eta_sd = 0.01)
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
  \item{pd0}{
%%     ~~Describe \code{pd0} here~~
}
  \item{i}{
%%     ~~Describe \code{i} here~~
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
function (eta.old, tree, tree.ht, pdist.old, pd0, i, Z, ywU, 
    eta_sd = 0.01) 
{
    eta.prop = eta.old + rnorm(length(eta.old), 0, sd = eta_sd)
    dist = 1/cophenetic(eb.phylo(tree, tree.ht, eta.prop))
    diag(dist) <- 0
    pdist.new = c(crossprod(dist[i, ], Z))
    no0 = which(pd0)
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
    ratio = min(1, exp(likeli))
    ratio
    u = (runif(1) <= ratio)
    if (u) {
        eta.old = eta.prop
        pdist.old = pdist.new
    }
    list(eta = eta.old, dist = pdist.old)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
