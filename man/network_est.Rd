\name{network_est}
\alias{network_est}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
network_est(Z, slices = 10, tree = NULL, model.type = c("full", "distance", "affinity"), uncertainty = FALSE, sparse = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Z}{
%%     ~~Describe \code{Z} here~~
}
  \item{slices}{
%%     ~~Describe \code{slices} here~~
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
  \item{sparse}{
%%     ~~Describe \code{sparse} here~~
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
function (Z, slices = 10, tree = NULL, model.type = c("full", 
    "distance", "affinity"), uncertainty = FALSE, sparse = TRUE, 
    ...) 
{
    require(geiger)
    require(phangorn)
    require(Matrix)
    require(methods)
    model.type = tolower(model.type[1])
    if (missing(Z)) 
        stop("Interaction matrix is missing!")
    if (slices == 0) 
        stop("no. of slices cannot be 0!")
    cleaned = network_clean(Z, tree, model.type)
    Z = cleaned$Z
    tree = cleaned$tree
    if (grepl("aff", model.type)) {
        print("Running affinity model...")
        param = fullJoint_est(Z, iter = slices, uncertainty = uncertainty, 
            ...)
    }
    if (grepl("(dist|full)", model.type)) {
        print(paste0("Running ", ifelse(grepl("dist", model.type), 
            "distance model...", "full model...")))
        param = ICM_est_over_acc(unname(Z), tree, slices, distOnly = grepl("dist", 
            model.type), uncertainty = uncertainty, sparse = sparse, 
            ...)
    }
    list(param = param, Z = Z, tree = tree)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
