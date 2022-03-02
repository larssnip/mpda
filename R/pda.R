#' @name pda
#' @title Pattern recognition with PLS+LDA
#'
#' @description A classification method that uses first PLS for dimension-reduction
#' and then LDA in the truncated score-space.
#'
#' @param y Vector of responses, must be a factor with exactly 2 levels. See \code{\link{mpda}}
#' for multi-level problems.
#' @param X Numeric matrix of predictor values.
#' @param prior Vector of prior probabilities, one value for each factor level in \code{y}.
#' @param max.dim Integer, the maximum number of dimensions to consider in PLS.
#' @param selected Vector of logicals, indicating a variable selection, see below.
#'
#' @details This classification method is designed for highly multivariate problems, i.e. where the
#' predictor matrix \code{X} has many and/or highly correlated columns (variables).
#'
#' First, the response factor is dummy-coded as 0's and 1's. This vector is then used together
#' with \code{X} to fit a PLS-model using the \code{oscorespls} algorithm, see the \code{\link{plsr}}
#' for details. The idea is that PLS will find linear combinations, denoted PLS-components, of the
#' original variables to be as an orthogonal basis for spanning the predictor space in such a way that objects
#' from the two factor levels are separated as much as possible. The score-matrix from this step are the original
#' data objects transformed into this subspace.
#'
#' Next, the score-matrix is truncated, i.e. only \code{max.dim} dimensions are used. The PLS-components
#' are all ordered such that the first component has the largest linear discriminative power. Thus, only a small
#' subspace is usually needed for separating between the two classes. This truncated score-matrix is used
#' as the predictor-matrix in LDA. One LDA-model is fitted for each dimension 1,...,\code{max.dim},
#' see \code{\link{lda}} for details.
#'
#' The predictor matrix \code{X} is centered, but not scaled. If you want scaled variables you need to do this
#' (with \code{\link{scale}}) before you call \code{pda}.
#'
#' The argument \code{selected} may be used to select a subset of the predictor variables (columns of \code{X}),
#' e.g. after a variable selection (see \code{\link{eliminator}}). This must be a
#' vector of logicals (\code{TRUE/FALSE}) indicating the selected variables, and the reduced predictor
#' matrix becomes \code{X[,which(selected)]}. The main reason for this option is the use of \code{\link{pda}}
#' in \code{\link{mpda}}.
#'
#' @return A \code{pda} object, which is a list with elements \code{PLS}, \code{LDA}, \code{Response} and
#' \code{Selected}. The
#' element \code{PLS} is simply the object returned from \code{\link{plsr}}. The element \code{LDA} is a
#' list with the fitted \code{lda} objects for each dimension. The elements \code{Response} and \code{Selected}
#' are copies of the arguments \code{y} and \code{selected}.
#'
#' @author Lars Snipen.
#'
#' @seealso \code{\link{predict.pda}}, \code{\link{pdaDim}}.
#'
#' @examples
#' data(microbiome)
#' y <- microbiome[1:40, 1]
#' X <- as.matrix(microbiome[1:40, -1])
#' m.trn <- pda(y, X, prior = c(0.5,0.5), max.dim = 10)
#'
#' data(poems)
#' y <- factor(poems[11:28,1], levels = c("Blake","Eliot"))
#' X <- as.matrix(poems[11:28, -1])
#' selection <- rep(FALSE, ncol(X))
#' selection[c(1,5,9,15,21)] <- TRUE   # using letters a, e, i, o and u only
#' p.trn <- pda(y, X, prior = c(1,1), selected = selection)
#'
#' @importFrom pls plsr
#' @importFrom MASS lda
#'
#' @export pda
#'
pda <- function(y, X, prior = NULL, max.dim = NULL, selected = NULL){
  if(!is.null(selected)){
    if(length(selected) != ncol(X)) stop("Argument selected must have ncol(X) elements")
    X <- X[, which( selected ), drop= F]
  }
  N <- nrow(X)
  P <- ncol(X)
  if(is.null(max.dim)) max.dim <- min(N-1, P-1)
  max.dim <- min(N, P, max.dim)
  if(!is.factor(y)) y <- factor(y)
  lev <- levels(y)
  L <- length(lev)
  if(L != 2) stop("Cannot have more than 2 factor levels in pda, use mpda for multi-level classification")
  y.dum <- as.numeric(y) - 1

  # The PLS step
  pls.mod <- plsr(y.dum ~ X, ncomp = max.dim, method = "oscorespls", validation = "none")
  Z <- pls.mod$scores

  # The LDA step
  if(is.null(prior)){
    prior <- as.numeric(table(y)/length(y))
  } else {
    prior <- prior/sum(prior)
  }
  lda.mod.lst <- vector("list", max.dim)
  for(i in 1:max.dim){
    lda.mod.lst[[i]] <- lda(Z[, 1:i, drop = F], y, prior = prior, tol = 1e-12)
  }
  pda.mod <- list(PLS = pls.mod, LDA = lda.mod.lst, Response = y, Selected = selected)
  class(pda.mod) <- c("pda", "list")
  return(pda.mod)
}




#' @name predict.pda
#' @title Classify based on a pda object
#'
#' @description Classify new data based on a trained \code{pda} model.
#'
#' @param object A fitted \code{\link{pda}} model.
#' @param newdata Matrix of predictor values.
#' @param ... Additional arguments to \code{\link{lda}}.
#'
#' @details Based on the trained \code{\link{pda}} model, new data objects (rows of \code{newdata}) are
#' classified. Remember that \code{\link{pda}} does not scale the predictor matrix, make certain your
#' \code{newdata} are treated identically to the predictor matrix used to train the \code{\link{pda}}
#' model.
#'
#' @return A list with one element for each dimension. In each element is another list with the integer
#' \code{Dimension}, the vector \code{Classifications}, the posterior probabilities \code{Posteriors} and
#' the PLS-scores \code{Scores} used to make the classifications.
#'
#' @author Lars Snipen.
#'
#' @seealso \code{\link{pda}}.
#'
#' @examples
#' data(microbiome)
#' y <- microbiome[1:40,1]
#' X <- as.matrix(microbiome[1:40, -1])
#' p.trn <- pda(y[-1], X[-1,], prior = c(0.5,0.5), max.dim = 8)    # leaving out sample 1 before training
#' lst <- predict(p.trn, newdata = X[1, , drop = FALSE])           # predicting sample 1
#'
#'
#' @importFrom stats predict
#'
#' @method predict pda
#' @export
predict.pda <- function(object, newdata=NULL, ... ){
  if(is.null(newdata)){
    newdata <- object$PLS$model$X
  } else {
    if(!is.null(object$Selected)){
      newdata <- newdata[, which(object$Selected), drop = F]
    }  # else use all columns
  }

  # The PLS step
  Z <- predict(object$PLS, newdata, type = "scores")
  p <- object$PLS$ncomp

  # The LDA step
  pda.lst <- vector("list", p)
  for(i in 1:p){
    lda.hat <- predict(object$LDA[[i]], Z[, 1:i, drop = F], ...)
    pda.lst[[i]] <- list(Dimension = i,
                         Classifications = lda.hat$class,
                         Posteriors = lda.hat$posterior,
                         Scores = Z[,1:i,drop = F])
  }
  return( pda.lst )
}




#' @name plot.pda
#' @title Plotting and summary of pda object
#' @aliases plot.pda summary.pda
#'
#' @description Generic functions for plotting and printing the content of a \code{pda} object.
#'
#' @param x A \code{pda} object, see below.
#' @param y not used.
#' @param object A \code{pda} object, see below.
#' @param col Two colors, one for each class-label.
#' @param pch Two markers, one for each class-label.
#' @param legend.pos Position of legend.
#' @param xlab Text on x-axis.
#' @param ylab Text on y-axis.
#' @param ... Optional graphical arguments.
#'
#' @details A \code{pda} object contains a fitted \code{\link{pda}} model.
#'
#' The \code{plot.pda} function display the samples as markers in the first 2 dimensions of the PLS-score space, and color
#' the markers by the class label information. If the PLS-score-space has only 1 dimension, the second axis in the plot
#' is non-informative. Since a \code{\link{pda}} object always has only 2 classes, you always specify a pair of colors
#' and markers in \code{col} and \code{pch}. The \code{legend.pos} must be one of the texts "bottomright", "bottom",
#' "bottomleft", "left", "topleft", "top", "topright", "right" and "center", specifying the position of the legend
#' inside the plot.
#'
#' The \code{summary.pda} function will display a text giving the number of samples from each class,
#' the number of PLS-dimensions used, and the priors used in the model fitting.
#'
#' @author Lars Snipen.
#'
#' @seealso \code{\link{pda}}.
#'
#' @examples
#' data(microbiome)
#' y <- microbiome[1:40,1]
#' X <- as.matrix(microbiome[1:40,-1])
#' p.trn <- pda(y,X,prior=c(0.5,0.5),max.dim=2)
#' plot(p.trn)
#' summary(p.trn)
#'
#' @importFrom graphics par plot legend
#'
#' @method plot pda
#' @export
plot.pda <- function(x, y = NULL, col = c("tan3","slategray3"), pch = c(15,15), legend.pos = "topright",
                     xlab = "PLS dimension 1", ylab = "PLS dimension 2", ...){
  Z <- x$PLS$scores
  n <- nrow(Z)
  p <- ncol(Z)
  y <- as.character(x$Response)
  lev <- levels(x$Response) #unique( y )
  cols <- col[as.numeric(x$Response)] #factor( y ) )]
  xx <- Z[,1]
  if(p > 1){
    yy <- Z[,2]
    y.axt = "s"
  } else {
    yy <- as.numeric(x$Response)
    ylab <- ""
    y.axt <- "n"
  }
  plot(xx, yy, pch = pch, col = cols, xlab = xlab, ylab = ylab, yaxt = y.axt, ...)
  legend(x = legend.pos, legend = lev, col = col, pch = pch)
}

#' @rdname plot.pda
#' @method summary pda
#' @export
summary.pda <- function(object, ...){
  lev <- levels(object$Response)
  cat("Fitted pda model with responses: ", lev[1], " (", sum(object$Response == lev[1]), ") and ",
      lev[2], " (", sum(object$Response == lev[2]), ")\n", sep = "")
  cat("using", length(object$LDA), "PLS-dimensions and priors:", object$LDA[[1]]$prior, "\n")
}








