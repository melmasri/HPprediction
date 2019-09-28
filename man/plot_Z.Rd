\name{plot_Z}
\alias{plot_Z}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot_Z(Z, xlab, ylab, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Z}{
%%     ~~Describe \code{Z} here~~
}
  \item{xlab}{
%%     ~~Describe \code{xlab} here~~
}
  \item{ylab}{
%%     ~~Describe \code{ylab} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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
function (Z, xlab, ylab, ...) 
{
    if (missing(ylab)) 
        ylab = "hosts"
    if (missing(xlab)) 
        xlab = "parasites"
    par(mar = c(5, 5, 1, 1) + 0.1)
    image(1:ncol(Z), 1:nrow(Z), t(Z[nrow(Z):1, ]), col = c("white", 
        "black"), ylab = ylab, xlab = xlab, useRaster = TRUE, 
        srt = 45, axes = FALSE, cex.lab = 3)
    a = 100 * max(1, round(ncol(Z) * 0.25/100, 0))
    axis(1, at = a * 0:(ceiling(ncol(Z)/a)), cex.axis = 2)
    b = 100 * max(1, round(nrow(Z) * 0.25/100, 0))
    axis(2, at = b * 0:(ceiling(nrow(Z)/b)), cex.axis = 2)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
