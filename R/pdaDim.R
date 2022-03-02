#' @name pdaDim
#' @title Finding optimal dimensionality in pda
#'
#' @description Fits \code{\link{pda}} models of 1,...,\code{max.dim} dimensions and finds the optimal
#' dimension by cross-validation and a regularization based on the McNemar-test.
#'
#' @param y Vector of responses, must be a factor with exactly 2 levels.
#' @param X Matrix of predictor values.
#' @param reg The regularization parameter, see below.
#' @param prior Vector of prior probabilities, one value for each factor level in \code{y}.
#' @param max.dim Integer, the maximum number of dimensions to consider.
#' @param selected Vector of logicals, indicating a variable selection, see below.
#' @param n.seg Integer, the number of cross-validation segments.
#' @param verbose Logical, turns on/off output during computations.
#'
#' @details The PLS method, which is part of the \code{\link{pda}} method, requires a decision on the
#' number of dimensions (components) to use. This is usually found by cross-validation, trying out many different
#' dimensions and searching for the best performance, e.g. the classification accuracy.
#'
#' In this algorithm a search for maximum accuracy is also conducted, but then the smallest dimension giving
#' an accuracy not significantly poorer than the maximum is used. The latter is based on the
#' McNemar-test, comparing the classifications of the maximum to that of the reduced model.
#' See \code{\link{mcnemar.test}} for details. This procedure is inspired by the CVANOVA idea of Indahl & Naes (1998).
#'
#' This implementation splits data into cross-validation segments, fits a \code{\link{pda}} model
#' and classifies. The accuracy for each dimension is the fraction of correctly classified elements
#' after the cross-validation.
#'
#' The McNemar-test is used to determine if a simpler model is not significantly poorer than the more complex
#' giving the maximum accuracy. The argument \code{reg} is the rejection level of this test, i.e. using
#' a \code{reg} value close to 0.0 means a harder regularization and a smaller dimension is in general selected.
#' Setting it to 1.0 means the dimension with the maximum accuracy (no regularization) is selected.
#'
#' The argument \code{selected} is used for variable selection, and just passed on to \code{\link{pda}}.
#'
#' @return A list with the elements \code{Dimension} and \code{Corrects}. \code{Dimension} is the
#' number of dimensions selected by this algorithm (integer). \code{Corrects} is a matrix of logicals
#' indicating which elements of \code{y} are correctly classified (\code{TRUE}) for each
#' dimension 1,...,\code{max.dim}.
#'
#' @author Lars Snipen.
#'
#' @references Indahl, UG, Naes, T (1998). Evaluation of alternative spectral feature
#' extraction methods of textural images for multivariate modeling. J. Chemometrics, 12:261-278.
#'
#' @seealso \code{\link{eliminator}}, \code{\link{mpda}}.
#'
#' @examples
#' data(microbiome)
#' y <- microbiome[1:40, 1]
#' X <- as.matrix(microbiome[1:40, -1])
#' lst <- pdaDim(y, X, reg = 0.1, prior = c(0.5,0.5), max.dim = 10)
#'
#' @importFrom pls plsr
#' @importFrom MASS lda qda
#' @importFrom stats mcnemar.test
#'
#' @export pdaDim
#'
pdaDim <- function(y, X, reg = 0.5, prior = NULL, max.dim = NULL, selected = NULL, n.seg = 10, verbose = TRUE){
  if(verbose) cat("pdaDim:\n")
  N <- nrow(X)
  P <- ncol(X)
  if(!is.null(selected)) P <- sum(selected)
  if(is.null(max.dim)) max.dim <- min(N, P)
  max.dim <- min(N, P, max.dim)
  y <- factor(y)
  ixx <- order(y)
  ys <- y[ixx]
  Xs <- X[ixx,,drop=F]
  seg <- rep(1:n.seg, length.out = N)

  # ensuring max.dim is small enough for cross-validation
  for(i in 1:n.seg){
    idx <- which(seg == i)
    ys.trn <- ys[-idx]
    Xs.trn <- Xs[-idx, , drop = F]
    md <- min(nrow(Xs.trn), ncol(Xs.trn))
    if(md < max.dim){
      max.dim <- md
      warning("max.dim is too large for cross-validation, set to ", max.dim, "\n")
    }
  }

  # the cross-validation
  if(verbose) cat("  cross-validation...\n")
  correct.mat <- matrix(rep(FALSE, N*max.dim), nrow = N)
  for(i in 1:n.seg){
    idx <- which(seg == i)
    ys.tst <- ys[idx]
    Xs.tst <- Xs[idx, , drop = F]
    ys.trn <- ys[-idx]
    Xs.trn <- Xs[-idx, , drop = F]
    px <- pda(ys.trn, Xs.trn, prior = prior, max.dim = max.dim, selected = selected)
    lst <- predict(px, Xs.tst)
    for(k in 1:length(lst)){
      correct.mat[idx,k] <- (lst[[k]]$Classifications == ys.tst)
    }
  }
  acc <- colMeans(correct.mat)
  idx.dim <- which(acc == max(acc))[1]
  opt.dim <- idx.dim
  if(verbose) cat("   maximum accuracy", format(acc[opt.dim], digits = 4), "at", opt.dim, "dimensions...\n")
  if(opt.dim > 1){
    for(j in seq((idx.dim-1), 1, -1)){
      contig.tab <- table(factor(correct.mat[,idx.dim], levels = c(FALSE, TRUE)),
                          factor(correct.mat[,j], levels = c(FALSE, TRUE)))
      tst <- mcnemar.test(contig.tab)
      if(tst$p.value > reg) opt.dim <- j    # continue with smaller dimension
    }
  }
  if(verbose) cat("   selected accuracy", format(acc[opt.dim], digits = 4), "at", opt.dim, "dimensions\n")
  pdim <- list(Dimension = opt.dim, Corrects = correct.mat)

  return(pdim)
}
