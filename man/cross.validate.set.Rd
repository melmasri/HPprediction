\name{cross.validate.set}
\alias{cross.validate.set}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
cross.validate.set(Z, rate = 0.2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Z}{
%%     ~~Describe \code{Z} here~~
}
  \item{rate}{
%%     ~~Describe \code{rate} here~~
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
function (Z, rate = 0.2) 
{
    if (max(range(Z)) > 1) 
        Z[Z > 0] <- 1
    Zo = Z
    pairs = which(Z == 1, arr.ind = T)
    n = floor(rate * nrow(pairs))
    sampled = rep(0, nrow(pairs))
    repeat {
        r = sample(which(sampled == 0), 1)
        if (sum(Z[pairs[r, 1], ]) > 1 & sum(Z[, pairs[r, 2]]) > 
            2) {
            Z[pairs[r, 1], pairs[r, 2]] <- 0
            sampled[r] <- 1
            n = n - 1
        }
        if (n == 0) 
            break
    }
    print(sprintf("Actual cross-validation rate is \%0.3f", sum(Zo - 
        Z)/sum(Zo)))
    Z
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
