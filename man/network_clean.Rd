\name{network_clean}
\alias{network_clean}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
network_clean(Z, tree = NULL, model.type = c("full", "distance", "affinity"), uncertainty = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Z}{
%%     ~~Describe \code{Z} here~~
}
  \item{tree}{
%%     ~~Describe \code{tree} here~~
}
  \item{model.type}{
%%     ~~Describe \code{model.type} here~~
}
  \item{uncertainty}{
%%     ~~Describe \code{uncertainty} here~~
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
function (Z, tree = NULL, model.type = c("full", "distance", 
    "affinity"), uncertainty = FALSE) 
{
    require(geiger)
    require(phangorn)
    require(Matrix)
    model.type = tolower(model.type[1])
    if (missing(Z)) 
        stop("Interaction matrix is missing!")
    if (!all(range(Z) == c(0, 1)) & !uncertainty) {
        warning("Z is converted to binary!", immediate. = TRUE, 
            call. = FALSE)
        Z = 1 * (Z > 0)
    }
    if (any(colSums(Z) == 0)) {
        stop("Z has empty columns, please remove!", immediate. = TRUE, 
            call. = FALSE)
        Z = Z[, which(colSums(Z) > 0)]
    }
    if (grepl("aff", model.type)) {
        if (!is.null(tree)) 
            warning("affinity model is chosen; ignoring tree!", 
                immediate. = TRUE, call. = FALSE)
    }
    if (grepl("(dist|full)", model.type)) {
        if (is.null(tree)) 
            stop("distance model is chosen, but tree is null!")
        if (!is.phylo(tree)) 
            stop("tree must be a phylogeny tree, see gieger!")
        if (!all(tree$tip.label \%in\% rownames(Z))) {
            warning("not all species in tree exist in Z; missing are removed from tree!", 
                immediate. = TRUE, call. = FALSE)
            tree = drop.tip(tree, tree$tip.label[!(tree$tip.label \%in\% 
                rownames(Z))])
        }
        if (!all(rownames(Z) \%in\% tree$tip.label)) {
            warning("not all row-species in Z exist tree; missing are removed from Z!", 
                immediate. = TRUE, call. = FALSE)
            Z = Z[rownames(Z) \%in\% tree$tip.label, ]
        }
        if (!all(sapply(sapply(rownames(Z), function(r) which(r == 
            tree$tip.label)), length) == 1)) 
            stop("some row-species in Z exist more than once tree!")
        aux = cophenetic(tree)
        if (!all(rownames(aux) == rownames(Z))) {
            row.order <- sapply(rownames(aux), function(r) which(r == 
                rownames(Z)))
            print("Ordering the rows of Z to match tree...")
            Z = Z[row.order, ]
        }
        if (max(range(tree$edge.length)) > 1) {
            print("normalizing tree edges by the maximum pairwise distance!")
            tree$edge.length = tree$edge.length/(max(aux)/2)
        }
    }
    list(Z = Z, tree = tree)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
