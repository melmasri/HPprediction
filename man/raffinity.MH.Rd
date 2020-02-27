\name{raffinity.MH}
\alias{raffinity.MH}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
raffinity.MH(old, z, Ud, sig = 0.1, hyper)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{old}{
%%     ~~Describe \code{old} here~~
}
  \item{z}{
%%     ~~Describe \code{z} here~~
}
  \item{Ud}{
%%     ~~Describe \code{Ud} here~~
}
  \item{sig}{
%%     ~~Describe \code{sig} here~~
}
  \item{hyper}{
%%     ~~Describe \code{hyper} here~~
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
function (old, z, Ud, sig = 0.1, hyper) 
{
    a = hyper[1]
    b = hyper[2]
    epsilon = rnorm(length(old), 0, sd = sig)
    new = old + -2 * (old + epsilon <= 0) * epsilon + epsilon
    likeli = (z + a - 1) * log(new/old) - (new - old) * (b + 
        Ud)
    ratio = pmin(1, exp(likeli))
    u = 1 * (runif(length(old)) <= ratio)
    u * new + (1 - u) * old
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
