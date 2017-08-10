#' @name microbiome
#' @docType data
#' @title Human Microbiome data
#'
#' @description Example data set from the \href{Human Microbiome Project}{http://hmpdacc.org/}.
#'
#' @usage
#' data(microbiome)
#'
#' @format A data set of 100 samples. The first column \code{Body.fluid} is a vector of class labels
#' for 5 different body-fluids/sites.
#' The remaining columns is a predictor matrix \code{X} with 1746 predictor variables.
#'
#' @details The vector \code{Body.fluid} is a factor with 5 levels, indicating the
#' body-fluid/site of each sample. The remaining columns can be used as a predictor matrix. Each row
#' corresponds to a taxonomic profile from the Human Microbiome Project. Each column contains
#' readcounts for some genus. Each sample (row) has been normalized to relative values (all rows
#' sum to 1.0)
#'
#' @author Lars Snipen.
#'
#' @seealso \code{\link{poems}}.
#'
#' @examples
#' data(microbiome)
#' y <- microbiome[,1]
#' X <- as.matrix(microbiome[,-1])
#'
NULL


#' @name poems
#' @docType data
#' @title Poetry data
#'
#' @description Example data set containing a response vector \code{Author} and a predictor
#' matrix with character frequencies for 28 english poems.
#'
#' @usage
#' data(poems)
#'
#' @format A data set of 28 samples. The first column \code{Author} is a vector of class labels
#' and the remaining columns is a predictor matrix \code{X} with 30 predictor variables.
#'
#' @details The vector \code{Author} is a factor with 3 levels, Blake, Eliot and Shakespeare, indicating the
#' author of each poem (sample). The remaining columns can be used as a predictor matrix. Each row
#' contains the relative character frequencies of 30 characters counted in each poem.
#'
#' @author Lars Snipen.
#'
#' @seealso \code{\link{poems}}.
#'
#' @examples
#' data(poems)
#' y <- poems[,1]
#' X <- as.matrix(poems[,-1])
#'
NULL

