\name{cross.validate.fold}
\alias{cross.validate.fold}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
cross.validate.fold(Z, n = 10, min.per.col = 1, missing.pattern = c("random", "prop.to.col.sums"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Z}{
%%     ~~Describe \code{Z} here~~
}
  \item{n}{
%%     ~~Describe \code{n} here~~
}
  \item{min.per.col}{
%%     ~~Describe \code{min.per.col} here~~
}
  \item{missing.pattern}{
%%     ~~Describe \code{missing.pattern} here~~
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
function (Z, n = 10, min.per.col = 1, missing.pattern = c("random", 
    "prop.to.col.sums")) 
{
    missing.pattern = tolower(missing.pattern[1])
    if (max(range(Z)) > 1) 
        Z[Z > 0] <- 1
    pairs = which(Z == 1, arr.ind = T)
    colnames(pairs) <- c("row", "col")
    if (length(which(colSums(Z) < min.per.col)) > 0) {
        aux = which(pairs[, "col"] \%in\% which(colSums(Z) < min.per.col))
        if (length(aux)) 
            pairs = pairs[-aux, ]
    }
    colm = pmax(colSums(Z) - min.per.col, 0)
    size = floor(sum(colm)/n)
    gr = rep(size, n)
    if (sum(colm)\%\%size != 0) 
        gr[n] = gr[n] + sum(colm)\%\%size
    group.colm = rep(1:n, times = gr)[sample.int(sum(colm), sum(colm))]
    pair.list = numeric(sum(colm))
    for (i in 1:sum(colm)) {
        a = which(colm > 0)
        if (missing.pattern == "random") 
            b = a[sample.int(length(a), 1)]
        else if (missing.pattern == "prop.to.col.sums") 
            b = a[sample.int(length(a), 1, prob = colm[a]/sum(colm[a]))]
        else stop("missing pattern has to be specified from selection!")
        colm[b] = colm[b] - 1
        pair.list[i] <- b
    }
    pair.list = tapply(pair.list, group.colm, identity)
    gr.list = list()
    bank = c()
    for (i in 1:n) {
        a = table(pair.list[[i]])
        gr.rows = unlist(sapply(1:length(a), function(r) {
            b = which(pairs[, "col"] == as.numeric(names(a[r])))
            b = setdiff(b, bank)
            b[sample.int(length(b), a[r])]
        }))
        bank = c(bank, gr.rows)
        gr.list[[i]] <- cbind(gr.rows, i)
    }
    aux = do.call("rbind", gr.list)
    pairs = cbind(pairs[aux[, 1], ], gr = aux[, 2])
    print(sprintf("Actual cross-validation rate is \%0.3f", table(pairs[, 
        "gr"])/sum(1 * (Z > 0))))
    pairs[order(pairs[, "gr"]), ]
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
