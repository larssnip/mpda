#' @name mpda
#' @title Multi-level pattern recognition with PLS+LDA
#'
#' @description A multi-level classification method that fits one \code{\link{pda}}-model to
#' every pair of factor levels.
#'
#' @param y Vector of responses, must be a factor with 3 or more levels.
#' @param X Numeric matrix of predictor values.
#' @param reg The regularization parameter used by \code{\link{pdaDim}}, see below.
#' @param prior Vector of prior probabilities, one value for each factor level in \code{y}.
#' @param selected Matrix of logicals, indicating selected predictor variables, see below.
#' @param max.dim Integer, the maximum number of dimensions to consider in PLS.
#' @param n.seg Integer, the number of cross-validation segments in \code{\link{pdaDim}}.
#' @param verbose Logical, turns on/off output during computations.
#'
#' @details A multi-level problem means the response \code{y} is a factor with at least 3 levels, i.e.
#' there are at least 3 distinct class labels (texts, integers, etc) in \code{y}. If you have a 2-level problem,
#' use \code{\link{pda}}.
#'
#' For each pair of factor levels, a \code{\link{pda}} model is fitted, using the subset of elements in \code{y}
#' containing these two levels, and the corresponding rows of \code{X} as predictors. If there are L factor levels,
#' there are N=L*(L-1)/2 such model-pairs.
#'
#' If the argument \code{reg} is positive (between 0 and 1) it means the \code{\link{pdaDim}} function is used
#' to estimate the dimensionality of each model. See \code{\link{pdaDim}} for details on this. If \code{reg} is
#' negative, all models will be fitted with \code{max.dim} dimensions.
#'
#' The argument \code{prior}, if specified, indicates the prior probability of each factor level. There must be one value for
#' each factor level, and in the exact same order as the factor levels. By default, all factor levels have equal prior.
#'
#' The matrix \code{selected} is used to indicate variable selection. For each pair of factor levels, a \code{\link{pda}}
#' model is fitted, and only a selection of the predictor variables (columns) in \code{X} may be used. In row k of
#' \code{selected} is a logical vector, one \code{TRUE/FALSE} value for each column in \code{X}. The \code{TRUE} elements
#' indicate which columns of \code{X} to use (\code{X[,which(selected[k,])]}). PLEASE NOTE: The rows of \code{selected}
#' must be ordered according to the factor level pairs:
#'
#' Row 1: Factor level 1 versus 2
#'
#' Row 2: Factor level 1 versus 3
#'
#' ...<all unique pairs of ordered levels>...
#'
#' Row N: Factor level (L-1) versus L
#'
#' @return An \code{mpda} object, which is a list of the pairwise \code{\link{pda}} objects. In addition,
#' the \code{mpda} object has two attributes: A copy of the full predictor matrix \code{X}
#' (\code{attr(mpda.obj,"X")})and the factor levels of \code{y} (\code{attr(mpda.obj,"Levels")}).
#'
#' @author Lars Snipen.
#'
#' @seealso \code{\link{predict.mpda}}, \code{\link{pdaDim}}.
#'
#' @examples
#' data(poems)
#' y <- poems[,1]
#' X <- as.matrix(poems[,-1])
#' mp.trn <- mpda(y,X,prior=c(1,1,1),max.dim=10)
#'
#' @importFrom pls plsr
#' @importFrom MASS lda
#'
#' @export mpda
#'
mpda <- function( y, X, reg=0.5, prior=NULL, selected=NULL, max.dim=NULL, n.seg=10, verbose=TRUE ){
  if( verbose ) cat( "mdpa: " )
  N <- nrow( X )
  P <- ncol( X )
  if( is.null( max.dim ) ) max.dim <- min( N, P )
  if( !is.factor(y) ) y <- factor( y )
  lev <- levels( y )
  L <- length( lev )
  N.pairs <- L*(L-1)/2
  if( L < 3 ) stop( "Cannot have less than 3 factor levels in mpda, use pda for 2-class classifications" )
  if( verbose ) cat( "Response with", L, "levels:", lev, "\n" )
  if( is.null( prior) ){
    prior <- as.numeric( table( y )/length( y ) )
  } else {
    if( length( prior ) != L ) stop( "Must have a prior value for each factor level" )
    prior <- prior/sum( prior )
    names( prior ) <- lev
    if( verbose ) cat( "The priors:", prior, "\n" )
  }
  if( !is.null( selected ) ){
    if( nrow( selected ) != N.pairs ) stop( "Matrix selected must have one row for each factor level pair in selected" )
    if( ncol( selected ) != P ) stop( "Matrix selected must have the same number of columns as X" )
  }

  cc <- 0
  selection <- NULL
  mpda.mod <- vector( "list", N.pairs )
  for( i in 1:(L-1) ){
    for( j in (i+1):L ){
      cc <- cc + 1
      if( verbose ) cat( "   fitting", lev[i], "versus", lev[j], "...\n" )
      idx <- which( y==lev[i] | y==lev[j] )
      yy <- factor( y[idx], levels=c(lev[i], lev[j]) )
      XX <- X[idx,]
      md <- min( nrow(XX), ncol(XX), max.dim )
      if( !is.null( selected ) ){
        selection <- selected[cc,]
        md <- min( md, sum( selection ) )
      }
      if( reg > 0 ){
        lst <- pdaDim( yy, XX, reg=reg, prior=prior[c(i,j)], max.dim=md, selected=selection, n.seg=n.seg, verbose=F )
        md <- lst$Dimension
      }
      pda.mod <- pda( yy, XX, prior=prior[c(i,j)], max.dim=md, selected=selection )
      mpda.mod[[cc]] <- pda.mod
    }
  }
  class( mpda.mod ) <- c( "mpda", "list" )
  attr( mpda.mod, "X" ) <- X
  attr( mpda.mod, "Levels" ) <- lev
  return( mpda.mod )
}



#' @name plot.mpda
#' @title Plotting and summary of mpda object
#' @aliases plot.mpda summary.mpda
#'
#' @description Generic functions for plotting and printing the content of a \code{mpda} object.
#'
#' @param x An \code{mpda} object, see below.
#' @param y not used.
#' @param object An \code{mpda} object, see below.
#' @param col Colors, one for each class-label.
#' @param pch Markers, one for each class-label.
#' @param legend.pos Position of legend.
#' @param xlab Text on x-axis.
#' @param ylab Text on y-axis.
#' @param ... Optional graphical arguments.
#'
#' @details An \code{mpda} object contains a fitted \code{\link{mpda}} model.
#'
#' The \code{plot.mpda} function produces a multi-panel plot of all the pairwise \code{\link{pda}} models,
#' see \code{\link{plot.pda}} for details
#'
#' The \code{summary.mpda} function will display a text giving the number of samples from each class,
#' the number of PLS-dimensions used, and the priors used in the model fitting.
#'
#' @author Lars Snipen.
#'
#' @seealso \code{\link{mpda}}.
#'
#' @examples
#' data(poems)
#' y <- poems[,1]
#' X <- as.matrix(poems[,-1])
#' mp.trn <- mpda(y,X,prior=c(1,1,1),max.dim=10)
#' plot(mp.trn)
#' summary(mp.trn)
#'
#' @method plot mpda
#' @export
plot.mpda <- function( x, y=NULL, col=c("tan3","slategray3"), pch=15, legend.pos="bottomright",
                      xlab="PLS dimension 1", ylab="PLS dimension 2", ... ){
  N <- length( x )
  p.rows <- floor( sqrt(N) )
  p.cols <- ceiling( N/p.rows )
  par( mfrow=c( p.rows, p.cols ) )
  for( i in 1:N ){
    plot( x[[i]], col=col, pch=pch, legend.pos=legend.pos, xlab=xlab, ylab=ylab )
  }
}

#' @rdname plot.mpda
#' @method summary mpda
#' @export
summary.mpda <- function( object, ... ){
  y <- unique( unlist( sapply( object, function(x){attr(x,"Classes")} ) ) )
  cat( "Fitted mpda model with responses:", y, "\n" )
}





#' @name predict.mpda
#' @title Classify based on an mpda object
#'
#' @description Classify new data based on a trained \code{mpda} model.
#'
#' @param object A fitted \code{\link{mpda}} model.
#' @param newdata Matrix of predictor values.
#' @param ... Additional arguments to \code{\link{lda}}.
#'
#' @details Based on the trained \code{\link{mpda}} model, new data objects (rows of \code{newdata}) are
#' classified, i.e. assigned to a level of the response factor. Remember that \code{\link{mpda}} does not scale
#' the predictor matrix, make certain your
#' \code{newdata} are treated identically to the predictor matrix used to train the \code{\link{mpda}}
#' model. If no \code{newdata} are available, the training data in the \code{mpda} object are used.
#'
#' Notice: Only the largest dimension of each pairwise \code{pda} model is used for prediction. When training
#' the \code{mpda} you would normally use the regularization argument \code{reg} and \code{\link{pdaDim}} to
#' estimate the dimension for each pariwise model, see \code{\link{mpda}}. Thus, different pairwise \code{pda}
#' are trained to their optimal dimensionality.
#'
#' @return A list with two elements, \code{Classifications} and \code{Post.means}.
#'
#' The vector \code{Classifications} are the suggested factor level for each sample in \code{newdata}.
#'
#' The matrix \code{Post.means} contains the posterior mean values for each factor level. For L factor levels,
#' each level is competing' against the other in (L-1) 'matches'. Each match results in posterior probabilities for the
#' two competing levels. The \code{Post.mean} scores of level k is simply the average of these (L-1) posterior probabilities.
#' If a sample has a \code{Post.means} value close to 1 for level k, it means this level is the clear winner in all
#' pairwise 'matches'.
#'
#' @author Lars Snipen.
#'
#' @seealso \code{\link{mpda}}.
#'
#' @examples
#' data(poems)
#' y <- poems[,1]
#' X <- as.matrix(poems[,-1])
#' mp.trn <- mpda(y,X,reg=0.5,prior=c(1,1,1),max.dim=3)
#' lst <- predict(mp.trn)
#' print(table(y,lst$Classifications))
#'
#' @importFrom stats predict
#'
#' @method predict mpda
#' @export
predict.mpda <- function( object, newdata=NULL, ... ){
  if( is.null( newdata ) ) newdata <- attr( object, "X" )
  N <- length( object )
  lev <- attr( object, "Levels" )
  L <- length( lev )
  Post.means <- matrix( 0, nrow=nrow( newdata ), ncol=L )
  colnames( Post.means ) <- lev
  rownames( Post.means ) <- rownames( newdata )
  for( i in 1:length( object ) ){
    lst <- predict( object[[i]], newdata=newdata )
    n <- length( lst )
    pm <- lst[[n]]$Posteriors
    idx <- which( lev == colnames( pm )[1] )
    Post.means[,idx] <- Post.means[,idx] + pm[,1]
    idx <- which( lev == colnames( pm )[2] )
    Post.means[,idx] <- Post.means[,idx] + pm[,2]
  }
  Post.means <- Post.means / (L-1)
  idx <- apply( Post.means, 1, function(x){which(x==max(x))[1]} )
  pred <- lev[idx]
  mpda.lst <- list( Classifications=pred, Post.means=Post.means )
  return( mpda.lst )
}
