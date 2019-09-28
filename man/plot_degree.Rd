\name{plot_degree}
\alias{plot_degree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot_degree(Z, Z_est, type = "both", host.col = "blue", parasite.col = "red")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Z}{
%%     ~~Describe \code{Z} here~~
}
  \item{Z_est}{
%%     ~~Describe \code{Z_est} here~~
}
  \item{type}{
%%     ~~Describe \code{type} here~~
}
  \item{host.col}{
%%     ~~Describe \code{host.col} here~~
}
  \item{parasite.col}{
%%     ~~Describe \code{parasite.col} here~~
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
function (Z, Z_est, type = "both", host.col = "blue", parasite.col = "red") 
{
    para_degrees <- as.data.frame(table(colSums(Z)))
    para_degrees$Var1 <- as.numeric(para_degrees$Var1)
    host_degrees <- as.data.frame(table(rowSums(Z)))
    host_degrees$Var1 <- as.numeric(host_degrees$Var1)
    xlim = c(1, max(para_degrees$Var1, host_degrees$Var1) * 1.5)
    ylim = c(1, max(para_degrees$Freq, host_degrees$Freq) * 1.5)
    if (!missing(Z_est)) {
        para_est <- as.data.frame(table(colSums(Z_est)))
        para_est$Var1 <- as.numeric(para_est$Var1)
        host_est <- as.data.frame(table(rowSums(Z_est)))
        host_est$Var1 <- as.numeric(host_est$Var1)
        xlim = c(1, max(para_degrees$Var1, host_degrees$Var1, 
            para_est$Var1, host_est$Var1) * 1.5)
        ylim = c(1, max(para_degrees$Freq, host_degrees$Freq, 
            para_est$Freq, host_est$Freq) * 1.5)
    }
    gpch = c("+", "*")
    if (type == "parasites") {
        plot((para_degrees), type = "p", col = parasite.col, 
            pch = gpch[2], log = "xy", xlim = xlim, ylim = ylim, 
            ylab = "Number of nodes", xlab = "Degree", cex.lab = 1.5, 
            cex.axis = 1.5)
        legend(xlim[2] * 0.03, ylim[2] * 1.4, c("Parasites"), 
            col = parasite.col, pch = gpch[2], bty = "n", pt.cex = 1.5, 
            cex = 2)
        if (!missing(Z_est)) {
            points((para_est), type = "p", col = parasite.col, 
                pch = 16)
            legend(xlim[2] * 0.03, ylim[2] * 0.8, c("Est"), col = parasite.col, 
                pch = 16, bty = "n", pt.cex = 1.5, cex = 2)
        }
    }
    if (type == "hosts") {
        plot((host_degrees), type = "p", col = host.col, pch = gpch[1], 
            log = "xy", xlim = xlim, ylim = ylim, ylab = "Number of nodes", 
            xlab = "Degree", cex.lab = 1.5, cex.axis = 1.5)
        legend(xlim[2] * 0.03, ylim[2] * 1.4, c("Hosts"), col = host.col, 
            pch = gpch[1], bty = "n", pt.cex = 1.5, cex = 2)
        if (!missing(Z_est)) {
            points((host_est), type = "p", col = host.col, pch = 16)
            legend(xlim[2] * 0.03, ylim[2] * 0.8, c("Est"), col = host.col, 
                pch = 16, bty = "n", pt.cex = 1.5, cex = 2)
        }
    }
    if (type == "both") {
        plot((para_degrees), type = "p", col = parasite.col, 
            pch = gpch[2], log = "xy", xlim = xlim, ylim = ylim, 
            ylab = "Number of nodes", xlab = "Degree", cex.lab = 1.5, 
            cex.axis = 1.5)
        points((host_degrees), type = "p", col = host.col, pch = gpch[1])
        legend(xlim[2] * 0.03, ylim[2] * 1.4, c("Parasites", 
            "Hosts"), col = c(parasite.col, host.col), pch = gpch[2:1], 
            bty = "n", pt.cex = 1.5, cex = 2)
        if (!missing(Z_est)) {
            points((para_est), type = "p", col = parasite.col, 
                pch = 16)
            points((host_est), type = "p", col = host.col, pch = 16)
            legend(xlim[2] * 0.03, ylim[2] * 0.33, c("Est"), 
                col = c("black"), pch = 16, bty = "n", pt.cex = 1.5, 
                cex = 2)
        }
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
