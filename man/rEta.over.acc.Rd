\name{rEta.over.acc}
\alias{rEta.over.acc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
rEta.over.acc(eta.old, pdist.old, pdist0, pdist00, sZ, Z, ywU, eta_sd = 0.01, sparse, dist.org, tmax)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eta.old}{
%%     ~~Describe \code{eta.old} here~~
}
  \item{pdist.old}{
%%     ~~Describe \code{pdist.old} here~~
}
  \item{pdist0}{
%%     ~~Describe \code{pdist0} here~~
}
  \item{pdist00}{
%%     ~~Describe \code{pdist00} here~~
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
function (eta.old, pdist.old, pdist0, pdist00, sZ, Z, ywU, eta_sd = 0.01, 
    sparse, dist.org, tmax) 
{
    eta.prop = eta.old + eta_sd * rnorm(1)
    dist = 1/EB.distance(dist.org, tmax, eta.prop)
    diag(dist) <- 0
    pdist.new = dist \%*\% sZ
    if (sparse) 
        pdist.new = as(pdist.new, "matrix")
    pdist.new[pdist00] <- 1
    acc = sapply(1:nrow(Z), function(i) {
        no0 = if (length(pdist0)) 
            pdist0[[i]]
        else NULL
        if (length(no0)) {
            if (length(no0) == sum(Z[i, ])) 
                likeli = -Inf
            else likeli = sum(((log(pdist.new[i, ]) - log(pdist.old[i, 
                ])) * Z[i, ])[-no0]) - sum((ywU[i, ] * (pdist.new[i, 
                ] - pdist.old[i, ]))[-no0])
        }
        else {
            likeli = sum((log(pdist.new[i, ]) - log(pdist.old[i, 
                ])) * Z[i, ]) - sum(ywU[i, ] * (pdist.new[i, 
                ] - pdist.old[i, ]))
        }
        likeli
    })
    acc = mean(pmin(1, exp(acc)))
    if (!is.nan(acc) && runif(1) <= acc) 
        return(list(eta = eta.prop, dist = pdist.new, change = TRUE))
    else return(list(eta = eta.old, change = FALSE))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
