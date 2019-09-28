\name{rExp2}
\alias{rExp2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
rExp2(l, g, Z, Z0, Z00)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{l}{
%%     ~~Describe \code{l} here~~
}
  \item{g}{
%%     ~~Describe \code{g} here~~
}
  \item{Z}{
%%     ~~Describe \code{Z} here~~
}
  \item{Z0}{
%%     ~~Describe \code{Z0} here~~
}
  \item{Z00}{
%%     ~~Describe \code{Z00} here~~
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
function (l, g, Z, Z0, Z00) 
{
    p = 1 - exp(-l)
    unif = matrix(runif(length(l)), dim(p))
    U = 1 + 0 * p
    U[-Z0] = -log(1 - unif[-Z0] * p[-Z0])/(l[-Z0] + tol.err)
    aa = Z00 & (unif < g * p/(tol.err + g * p + 1 - p))
    U[aa] = -log(1 - unif[aa] * (g * p[aa] + 1 - p[aa])/g)/l[aa]
    U
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
