\name{topPairs}
\alias{topPairs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
topPairs(P, Z, topX = 20)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{P}{
%%     ~~Describe \code{P} here~~
}
  \item{Z}{
%%     ~~Describe \code{Z} here~~
}
  \item{topX}{
%%     ~~Describe \code{topX} here~~
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
function (P, Z, topX = 20) 
{
    require(reshape2)
    P[Z > 0] <- -1
    aux = melt(P)
    aux = aux[order(aux$value, decreasing = TRUE), ]
    colnames(aux) <- c("Hosts", "Parasites", "p")
    aux[1:topX, ]
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
