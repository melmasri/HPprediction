\name{rEta.copheneticSuperFast}
\alias{rEta.copheneticSuperFast}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
rEta.copheneticSuperFast(eta.old, pdist.old, no0, i, sZ, Z, ywU, eta_sd = 0.01, sparse, dist.org, tmax)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eta.old}{
%%     ~~Describe \code{eta.old} here~~
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
  \item{sparse}{
%%     ~~Describe \code{sparse} here~~
}
  \item{dist.org}{
%%     ~~Describe \code{dist.org} here~~
}
  \item{tmax}{
%%     ~~Describe \code{tmax} here~~
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
function (eta.old, pdist.old, no0, i, sZ, Z, ywU, eta_sd = 0.01, 
    sparse, dist.org, tmax) 
{
    change = FALSE
    if (length(no0) && length(no0) == sum(Z[i, ])) {
        likeli = -Inf
    }
    else {
        eta.prop = eta.old + eta_sd * rnorm(1)
        dist = 1/EB.distance(dist.org[i, ], tmax, eta.prop)
        dist[i] <- 0
        pdist.new = dist \%*\% sZ
        if (sparse) 
            pdist.new = pdist.new@x
        if (length(no0)) {
            likeli = sum(((log(pdist.new) - log(pdist.old)) * 
                Z[i, ])[-no0]) - sum((ywU * (pdist.new - pdist.old))[-no0])
        }
        else {
            likeli = sum((log(pdist.new) - log(pdist.old)) * 
                Z[i, ]) - sum(ywU * (pdist.new - pdist.old))
        }
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
