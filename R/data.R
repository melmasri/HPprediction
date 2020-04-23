#' Host-Parasite interaction data from the Global Mammal Parasite Database 2.0
#'
#' The GMPD 2.0 is a database of the parasites of wild ungulates (artiodactyls and perissodactyls), 
#' carnivores, and primates. It comprises 24,000 entries representing data from over 2700 literature sources. 
#' Host Latin binomials (column "HostCorrectedName") match those used by the Fritz et al. (2009) mammal supertree (data(mammal_supertree)).
#' 
#' @docType data
#'
#' @usage data(gmpd)
#'
#' @format An object of class \code{"data.frame"}.
#'
#' @keywords datasets
#' 
#' @references Stephens, P.R., Pappalardo, P., Huang, S., Byers, J.E., 
#' Farrell, M.J., Gehman, A., Ghai, R.R., Haas, S.E., Han, B., Park, A.W., 
#' Schmidt, J.P., Altizer, S., Ezenwa, V.O. and Nunn, C.L. (2017), 
#' Global Mammal Parasite Database version 2.0. 
#' Ecology, 98: 1476-1476. \href{https://doi.org/10.1002/ecy.1799}{doi:10.1002/ecy.1799}
#'
#' @source \href{https://doi.org/10.1002/ecy.1799}{Ecology Data Papers}
#'
#' @examples
#' data(gmpd)
"gmpd"


#' Fritz et al. (2009) mammal supertree
#'
#' A phylogeny of 5020 extant mammals, with names harmonized to Wilson & Reeder (2005).
#' 
#' @docType data
#'
#' @usage data(mammal_supertree)
#'
#' @format An object of class \code{"phylo"}.
#'
#' @keywords datasets
#' 
#' @references Fritz, S.A., Bininda‐Emonds, O.R.P. and Purvis, A. (2009), 
#' Geographical variation in predictors of mammalian extinction risk: big is bad, but only in the tropics. 
#' Ecology Letters, 12: 538-549. \href{https://doi.org/10.1111/j.1461-0248.2009.01307.x}{doi:10.1111/j.1461-0248.2009.01307.x}
#' 
#' Wilson, D.E. & Reeder, D.M. 2005. Mammal Species of the World: A Taxonomic and Geographic Reference. Third Edition. The Johns Hopkins University Press, Baltimore..
#'
#' @source \href{https://doi.org/10.1111/j.1461-0248.2009.01307.x}{Ecology Letters}
#'
#' @examples
#' data(mammal_supertree)
"mammal_supertree"