\name{fullJoint_est}
\alias{fullJoint_est}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
fullJoint_est(Z, iter = 10, uncertainty = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Z}{
%%     ~~Describe \code{Z} here~~
}
  \item{iter}{
%%     ~~Describe \code{iter} here~~
}
  \item{uncertainty}{
%%     ~~Describe \code{uncertainty} here~~
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
function (Z, iter = 10, uncertainty = FALSE, ...) 
{
    el <- list(...)
    nw = ncol(Z)
    ny = nrow(Z)
    n = nw * ny
    y = if (!is.null(el$y)) 
        el$y
    else 1
    w = if (!is.null(el$w)) 
        el$w
    else 1
    a_w = if (!is.null(el$a_w)) 
        el$a_w
    else 1
    a_y = if (!is.null(el$a_y)) 
        el$a_y
    else 1
    b_w = if (!is.null(el$b_w)) 
        el$b_w
    else 1
    b_y = if (!is.null(el$b_y)) 
        el$b_y
    else 1
    y_sd = if (is.null(el$y_sd)) 
        0.2
    else el$y_sd
    w_sd = if (is.null(el$w_sd)) 
        0.2
    else el$w_sd
    batch.size = if (!is.null(el$batch.size)) 
        el$batch.size
    else 50
    beta = if (!is.null(el$beta)) 
        el$beta
    else 1
    burn.in = if (is.null(el$burn.in)) 
        floor(0.5 * iter)
    else floor(el$burn.in * iter)
    y0 <- matrix(y, nrow = ny, ncol = iter + 1)
    w0 <- matrix(w, nrow = nw, ncol = iter + 1)
    g0 <- rep(0, iter + 1)
    peta = rep(0, iter + 1)
    Z0 <- Z == 0
    mc <- colSums(Z)
    mr <- rowSums(Z)
    print(sprintf("Run for \%i iterations with \%i burn-ins", iter, 
        burn.in))
    s = 1
    tryCatch(for (s in 1:iter) {
        if (s\%\%200 == 0) {
            print(sprintf("iteration \%d, at \%s", s, Sys.time()))
            if (!is.null(el$backup)) 
                save(y0, w0, peta, g0, file = "snapshot.RData")
        }
        U0 <- rExp(outer(y0[, s], w0[, s]))
        U0[Z0] <- 1
        w0[, s + 1] <- raffinity.MH(w0[, s], mc, crossprod(y0[, 
            s], U0), sig = w_sd, c(a_w, b_w))
        y0[, s + 1] <- raffinity.MH(y0[, s], mr, tcrossprod(w0[, 
            s + 1], U0), sig = y_sd, c(a_y, b_y))
        if (s\%\%batch.size == 0) {
            ss = s + 1 - 1:batch.size
            ac = 1 - rowMeans(abs(w0[, ss] - w0[, ss + 1]) < 
                tol.err)
            w_sd = w_sd * exp(beta * (ac - 0.44)/log(1 + s/batch.size))
            ac = 1 - rowMeans(abs(y0[, ss] - y0[, ss + 1]) < 
                tol.err)
            y_sd = y_sd * exp(beta * (ac - 0.44)/log(1 + s/batch.size))
        }
    }, warning = function(war) {
        print(c("warning at iter:", s))
        print(war)
        traceback()
    }, error = function(e) {
        print(c("error at iter:", s))
        print(e)
        traceback()
    }, finally = print("Done!"))
    if (burn.in == 0) 
        burn.in = 1
    else burn.in = 1:(burn.in + 1)
    y0 = y0[, -burn.in]
    w0 = w0[, -burn.in]
    if (uncertainty) 
        g0 = g0[-burn.in]
    else g0 = NULL
    list(w = w0, y = y0, g = g0, burn.in = max(burn.in) - 1, 
        sd = list(w = w_sd, y = y_sd))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
