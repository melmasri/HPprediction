\name{sample_parameter}
\alias{sample_parameter}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sample_parameter(param, MODEL, Z, tree, size = 1000, weights = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{param}{
%%     ~~Describe \code{param} here~~
}
  \item{MODEL}{
%%     ~~Describe \code{MODEL} here~~
}
  \item{Z}{
%%     ~~Describe \code{Z} here~~
}
  \item{tree}{
%%     ~~Describe \code{tree} here~~
}
  \item{size}{
%%     ~~Describe \code{size} here~~
}
  \item{weights}{
%%     ~~Describe \code{weights} here~~
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
function (param, MODEL, Z, tree, size = 1000, weights = NULL) 
{
    sample_mcmc <- function(mcmc_sample, nObs, size = 1000) {
        if (is.matrix(mcmc_sample)) 
            mcmc_sample[, sample.int(nObs, size, replace = TRUE)]
        else mcmc_sample[sample.int(nObs, size, replace = TRUE)]
    }
    if (!is.null(weights) & length(weights) != size) 
        stop("weights are not the same sampling size.")
    t.max = get.max.depth(tree)
    dist.original = unname(cophenetic(rescale(tree, "EB", 0)))/2
    if (grepl("dist", MODEL)) {
        Y = W = 1
    }
    else {
        Y = sample_mcmc(param$y, ncol(param$y), size)
        W = sample_mcmc(param$w, ncol(param$w), size)
    }
    if (grepl("(dist|full)", MODEL)) {
        Eta = sample_mcmc(param$eta, length(param$eta), size)
    }
    else Eta = 1
    zeroZ = which(Z > 0)
    P <- 0
    for (s in 1:size) {
        if (grepl("(aff|full)", MODEL)) {
            YW = outer(Y[, s], W[, s])
        }
        else YW = 1
        if (grepl("(full|dist)", MODEL)) {
            distance = 1/EB.distance(dist.original, t.max, Eta[s])
            diag(distance) <- 0
            distance = distance \%*\% Z
            distance[distance == 0] <- if (grepl("dist", MODEL)) 
                Inf
            else 1
        }
        else distance = 1
        if (!is.null(weights)) {
            a = 1 - exp(-YW * distance)
            Pg = a * weights[s]/(1 - a + weights[s] * a)
            Pg[zeroZ] <- a[zeroZ]
            P = P + Pg
        }
        else {
            P = P + 1 - exp(-YW * distance)
        }
    }
    matrix(P/size, nrow = nrow(com), ncol = ncol(com))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
