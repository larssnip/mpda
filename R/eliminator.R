#' @name eliminator
#' @title The Eliminator algorithm
#'
#' @description Variable elimination using pda.
#'
#' @param y Vector of responses, a factor of exact 2 levels.
#' @param X Matrix of predictor values.
#' @param reg The regularization parameter, see below.
#' @param prior Vector of prior probabilities, one value for each factor level in \code{y}.
#' @param max.dim Integer, the maximum number of dimensions to consider.
#' @param frac Fraction of unimportant variables to eliminate in each iteration (default is 0.25).
#' @param vip.lim The threshold for the VIP criterion (default is 1.0).
#' @param n.seg Integer, the number of cross-validation segments (default is 10).
#' @param verbose Logical, turns on/off output during computations.
#'
#' @details This is a variable selection (elimination) algorithm based on the \code{\link{pda}} method
#' where the response (\code{y}) is a factor with 2 levels, i.e. a two-class problem. The restriction
#' to a two-class problem comes from the use of the VIP criterion. The idea is a slight modification
#' of the one described in Mehmood et al, (2011).
#'
#' The algorithm starts out by using all variables, fitting a \code{\link{pda}} model, using the
#' \code{\link{pdaDim}} function to estimate the optimal dimensionality for the PLS step after each
#' elimination. The regularization parameter \code{reg} and the number of
#' cross-validation segments \code{n.seg} are needed for this, see \code{\link{pdaDim}} for details.
#'
#' Based on the fitted model, the VIP-criterion is computed for all variables, see \code{\link{vipCriterion}}
#' for details. Variables are ranked by this criterion, and those with a VIP-score below \code{vip.lim}
#' are considered unimportant. The fraction \code{frac} of the unimportant variables are eliminated.
#' If no unimportant variables are left, the algorithm terminates. It will mostly continue until
#' only 1 variable is left, and then terminate.
#'
#' After each elimination, the classification accuracy, i.e. the fraction of correctly classified samples, is
#' computed based on the cross-validation of the \code{\link{pdaDim}} results. The reduced model accuracy is
#' tested against the maximum so far, using the \code{\link{mcnemar.test}}. The p-value of this test is
#' reported for each elimination step, indicating if
#' the reduced model has a significantly poorer accuracy than the maximum so far. The idea is to keep
#' on eliminating variables until the accuracy becomes significantly poorer than the maximum.
#'
#' The argument \code{max.dim} allows you to specify the maximum number of PLS components to consider. The
#' optimal dimension found by \code{\link{pdaDim}} is between 1 and \code{max.dim}. If it turns out identical
#' to your \code{max.dim}, increase its value slightly and re-run.
#'
#' @return A \code{list} with two matrices, \code{Elimination} and \code{Selected}.
#'
#' \code{Elimination} has one row for each iteration, with accuracy results after each. The first column is
#' the number of variables left, the second is the fraction of correctly classified samples, and the third
#' is the p-value of the McNemar-test, see above.
#' You will typically look down this matrix for iterations where you have the maximum accuracy. Then,
#' continue down the matrix from that row, until there is a significant drop in accuracy (column 2) and a
#' corresponding small p-value (column 3), and choose the result just
#' before this as your final selection.
#'
#' In \code{Selected} you also find one row for each iteration. This logical indicates which variables are
#' selected after each iteration.
#' Thus, if cell \code{[i,j]} is \code{TRUE} variable \code{j} is still included after iteration \code{i}.
#'
#' @author Lars Snipen.
#'
#' @references Mehmood, T, Martens, H, Saeboe, S, Warringer, J, Snipen, L (2011). A Partial Least Squares based algorithm
#' for parsimonious variable selection. Algorithms for Molecular Biology, 6:27.
#'
#' @seealso \code{\link{vipCriterion}}, \code{\link{pdaDim}}.
#'
#' @examples
#' data(microbiome)
#' y <- microbiome[1:40,1]
#' X <- as.matrix(microbiome[1:40,-1])
#' lst <- eliminator(y,X,max.dim=10)
#' print(lst$Elimination)
#' # Seems like iteration 23 is the place to stop, since a significant drop
#' # in performance is found at iteration 24.
#' # There are 2 variables left in iteration 23. These are variables
#' print(which(lst$Selected[23,]))
#'
#' @importFrom stats mcnemar.test
#'
#' @export eliminator
#'
eliminator <- function( y, X, reg=0.5, prior=NULL, max.dim=NULL,
                        frac=0.25, vip.lim=1.0, n.seg=10, verbose=TRUE ){
  if( verbose ) cat( "The Eliminator:\n" )
  N <- nrow( X )
  if( n.seg > N ) stop( "Must use n.seg less than nrow(X)!" )
  P <- ncol( X )
  if( P==1 ) stop( "Meaningless to perform variable-selection with 1 variable!" )
  y <- factor( y )
  if( nlevels( y ) > 2 ) stop( "Response y must be a factor of exactly 2 levels" )
  y <- as.integer( y )

  # All variables
  selected <- 1:P
  pdim <- pdaDim( y, X, reg=reg, prior=prior, max.dim=max.dim, n.seg=n.seg, verbose=F )
  c.max <- pdim$Corrects[,pdim$Dimension]
  if( verbose ) cat( "   full model has", P, "variables, accuracy =",
                     format( sum(c.max)/length(c.max), digits=4 ), "using", pdim$Dimension, "dimensions\n" )
  px <- pda( y, X, max.dim=pdim$Dimension )
  vips <- vipCriterion( px$PLS, pdim$Dimension )
  idx.unimportant <- which( vips < vip.lim )
  eliminate <- (length( idx.unimportant ) > 0) # TRUE if there are variables with VIP below vip.lim
  N.vars <- length( selected )
  correct <- sum( c.max )
  pvalues <- 1
  select.mat <- matrix( T, nrow=1, ncol=P )

  # The elimination
  zerow <- matrix( FALSE, nrow=1, ncol=P, dimnames=NULL )
  while( eliminate ){
    # Eliminating
    n <- ceiling( frac*length( idx.unimportant ) )
    ixx <- order( vips )
    selected <- selected[-ixx[1:n]]
    XX <- X[,selected,drop=F]
    PP <- ncol( XX )

    # Computing accuracy based on reduced variable set
    pdim <- pdaDim( y, XX, reg=reg, prior=prior, max.dim=max.dim, n.seg=n.seg, verbose=F )
    crrct <- pdim$Corrects[,pdim$Dimension]
    if( sum( crrct ) >= sum( c.max ) ){
      c.max <- crrct
      pvl <- 1
    } else {
      ct <- table( factor( c.max, levels=c(FALSE,TRUE) ), factor( crrct, levels=c(FALSE,TRUE) ) )
      tst <- mcnemar.test( ct )
      pvl <- tst$p.value
    }
    N.vars <- c( N.vars, PP )
    correct <- c( correct, sum( crrct ) )
    pvalues <- c( pvalues, pvl )
    newrow <- zerow
    newrow[selected] <- TRUE
    select.mat <- rbind( select.mat, newrow )
    if( verbose ) cat( "   eliminated to", PP, "variables, accuracy =",
                         format( sum(crrct)/length(crrct), digits=4 ), "using",
                         pdim$Dimension, "dimensions\n" )

    # Computing VIPs for the reduced model, for next step
    px <- pda( y, XX, max.dim=pdim$Dimension )
    vips <- vipCriterion( px$PLS, pdim$Dimension )
    idx.unimportant <- which( vips < vip.lim )
    eliminate <- (length( idx.unimportant ) > 0) & (PP > 1)
  }
  colnames( select.mat ) <- colnames( X )
  rownames( select.mat ) <- paste( "Iteration", 1:nrow( select.mat ) )
  emat <- matrix( c( N.vars, correct/length(y), pvalues ), ncol=3, byrow=F )
  colnames( emat ) <- c( "N.variables", "Accuracy", "P.value" )
  rownames( emat ) <- rownames( select.mat )

  return( list( Elimination=emat, Selected=select.mat ) )
}



#' @name mEliminator
#' @title The Eliminator for mpda
#'
#' @description Variable elimination in mpda.
#'
#' @param y Vector of responses, a factor of exact 2 levels.
#' @param X Matrix of predictor values.
#' @param reg1 The regularization parameter for \code{\link{pdaDim}}.
#' @param reg2 The regularization parameter for selection, see below.
#' @param prior Vector of prior probabilities, one value for each factor level in \code{y}.
#' @param max.dim Integer, the maximum number of dimensions to consider.
#' @param frac Fraction of unimportant variables to eliminate in each iteration (default is 0.25).
#' @param vip.lim The threshold for the VIP criterion (default is 1.0).
#' @param n.seg Integer, the number of cross-validation segments (default is 10).
#' @param verbose Logical, turns on/off output during computations.
#'
#' @details This is a wrapper for doing variable selection with the \code{\link{eliminator}} on an
#' \code{\link{mpda}} object.
#'
#' You use this function if you have a multi-level classification problem, and wants
#' a standardized (and regularized) variable selection. This function uses \code{\link{mpda}} for the
#' multi-level problem, which means all pairs of levels are modelled. A variable selection is performed
#' for each level-pair, using the \code{\link{eliminator}} algorithm.
#'
#' The argument \code{reg2} is a regularization parameter along the same line as \code{reg1}, which is used
#' by \code{\link{pdaDim}}. It is a rejection level of the \code{\link{mcnemar.test}}. In the \code{\link{eliminator}}
#' algorithm, this test is performed after each elimination step, to see if the resulting accuracy is significantly
#' pooerer than the maximum accuracy seen up to that step. As long as the corresponding p-value is at least as
#' large as \code{reg2}, the elimination should continue. Thus, setting \code{reg2=1.0} (default) means there is
#' no regularization, and the selection producing the maximum accuracy is the result. By lowering \code{reg2} you
#' get a more stable selection, at the potential cost of elimination too much.
#'
#' @return A matrix with one row for each level-pair and one column for each variable (column) in \code{X}.
#'
#' Each row is a logical vector indicating which variables (\code{TRUE}) that were selected for the corresponding
#' level-pair. Thus, if we denote this matrix \code{S}, then \code{X[,S[1,]]} is the sub-matrix of \code{X} selected
#' to be optimal for the use for level-pair \code{1}, etc.
#'
#' @author Lars Snipen.
#'
#' @seealso \code{\link{eliminator}}, \code{\link{mpda}}.
#'
#' @examples
#' data(poems)
#' y <- poems[,1]
#' X <- as.matrix(poems[,-1])
#' # Variable selection
#' S <- mEliminator(y,X,max.dim=10)
#'
#' # Fitting model with selection information
#' mp.trn <- mpda(y,X,prior=c(1,1,1),selected=S,max.dim=10)
#' # Predicting...
#' predict(mp.trn)
#'
#' @export mEliminator
#'
mEliminator <- function( y, X, reg1=0.5, reg2=1.0, prior=NULL, max.dim=NULL,
                         frac=0.25, vip.lim=1.0, n.seg=10, verbose=TRUE ){
  if( verbose ) cat( "mEliminator: " )
  N <- nrow( X )
  P <- ncol( X )
  if( !is.factor(y) ) y <- factor( y )
  lev <- levels( y )
  L <- length( lev )
  N.pairs <- L*(L-1)/2
  if( L < 3 ) stop( "mEliminator is made for multi-level problems, use the eliminator for 2-level problems" )
  if( verbose ) cat( "Response with", L, "levels:", lev, "\n" )
  if( is.null( prior) ){
    prior <- as.numeric( table( y )/length( y ) )
  } else {
    if( length( prior ) != L ) stop( "Must have a prior value for each factor level" )
    prior <- prior/sum( prior )
    names( prior ) <- lev
    if( verbose ) cat( "The priors:", prior, "\n" )
  }

  cc <- 0
  Selected <- matrix( TRUE, nrow=N.pairs, ncol=P )
  for( i in 1:(L-1) ){
    for( j in (i+1):L ){
      cc <- cc + 1
      if( verbose ) cat( "Eliminating for", lev[i], "versus", lev[j], "...\n" )
      idx <- which( y==lev[i] | y==lev[j] )
      yy <- factor( y[idx], levels=c(lev[i], lev[j]) )
      XX <- X[idx,]
      md <- min( nrow(XX), ncol(XX), max.dim )
      lst <- eliminator( yy, XX, reg=reg1, prior=prior[c(i,j)], max.dim=md,
                         frac=frac, vip.lim=vip.lim, n.seg=n.seg, verbose=verbose )
      idx.max <- max( which( lst$Elimination[,2] == max( lst$Elimination[,2] ) ) )
      idx.brk <- max( which( lst$Elimination[idx.max:nrow( lst$Elimination ),3] >= reg2 ) )
      idx.sel <- idx.max + idx.brk - 1
      Selected[cc,] <- lst$Selected[idx.sel,]
      if( verbose ) cat( "Selected", sum( Selected[cc,] ), "variables\n" )
    }
  }
  return( Selected )
}






