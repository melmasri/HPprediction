\name{rocCurves}
\alias{rocCurves}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
rocCurves(Z.test, Z.train, P, plot = TRUE, bins = 400, all = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Z.test}{
%%     ~~Describe \code{Z.test} here~~
}
  \item{Z.train}{
%%     ~~Describe \code{Z.train} here~~
}
  \item{P}{
%%     ~~Describe \code{P} here~~
}
  \item{plot}{
%%     ~~Describe \code{plot} here~~
}
  \item{bins}{
%%     ~~Describe \code{bins} here~~
}
  \item{all}{
%%     ~~Describe \code{all} here~~
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
function (Z.test, Z.train, P, plot = TRUE, bins = 400, all = FALSE) 
{
    Z.test = 1 * (Z.test > 0)
    Z.train = 1 * (Z.train > 0)
    u = seq(0, 1, length.out = bins)
    Z1 = (Z.test - Z.train) == 1
    if (all) 
        Z1 = Z.test == 1
    m = sum(1 * (Z1))
    Z2 = Z.test == 0
    n = sum(1 * Z2)
    aux = sapply(u, function(r) {
        aux = 1 * (P >= r)
        TP = sum(aux[Z1])
        FN = m - TP
        FP = sum(aux[Z2])
        TN = n - FP
        c(TP = TP, FN = FN, FP = FP, TN = TN)
    })
    FPR = aux["FP", ]/(aux["FP", ] + aux["TN", ])
    TPR = aux["TP", ]/(aux["TP", ] + aux["FN", ])
    roc = data.frame(u = u, FPR = FPR, TPR = TPR)
    max.point = which.min(abs(roc$TPR - 1) + roc$FPR)
    threshold = roc$u[max.point]
    n = length(u)
    auc = 0.5 * t(abs(FPR[2:n] - FPR[2:n - 1])) \%*\% (TPR[2:n] + 
        TPR[2:n - 1])
    auc = round(100 * auc, 2)
    if (plot) {
        plot(FPR, TPR, xlab = "FPR (1-specifity)", ylab = "TPR (sensitivity)", 
            type = "l", col = "red", lwd = 2, main = paste("AUC ", 
                auc), xlim = c(0, 1), ylim = c(0, 1), pch = 6, 
            lty = 4)
        abline(a = 0, b = 1, col = "black", lty = 2, lwd = 2)
    }
    list(auc = auc, threshold = threshold, roc = roc)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
