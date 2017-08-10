#' @name vipCriterion
#' @title The VIP criterion
#'
#' @description Computes the VIP criterion used to rank variable importance.
#'
#' @param pls.model Object with fitted PLS-model.
#' @param dim Integer, the number of dimensions to consider.
#'
#' @details After the fitting of a PLS-model, some of the original variables will have more
#' impact than the others on the prediction of the response. The VIP criterion is one way
#' to quantify this, see Chong&Yun, 2005. This criterion requires a single response regression
#' problem, which means a two-class classification problem.
#'
#' A large VIP indicates the corresponding variable is important. A threshold at 1.0 is often
#' used, variables with VIP above 1.0 are the important ones.
#'
#' @return A vector of VIP scores, one for each variable in the predictor matrix of the fitted
#' PLS model.
#'
#' @author Lars Snipen.
#'
#' @references Chong, Il-Gyo & Jun, Chi-Hyuck, 2005, Performance of some variable selection
#' methods when multicollinearity is present, Chemometrics and Intelligent Laboratory
#' Systems 78, 103--112.
#'
#' @seealso \code{\link{eliminator}}.
#'
#' @export vipCriterion
#'
vipCriterion <- function( pls.model, dim=1 ){
  if( nrow(pls.model$Yloadings) > 1 ) stop( "Only implemented for single-response models" )
  W <- pls.model$loading.weights[,1:dim,drop=F]    # (P x a) matrix
  Z <- pls.model$scores[,1:dim,drop=F]             # (N x a) matrix
  q  <- as.numeric( pls.model$Yloadings )[1:dim]   # (1 x a) vector
  P  <- nrow( W )                                  # scalar

  explained.variance <- q^2 * colSums( Z^2 )
  W.norm <- colSums( W^2 )     # Usually not needed, this norm is 1.0 anyway (for most PLS algorithms)
  importance <- W^2 / W.norm   # Importance = squared loading weights
  denominator <- sum( explained.variance )
  numerator <- rowSums( matrix( rep( explained.variance, P ), nrow=P, byrow=T ) * importance )

  vip <- sqrt( P * numerator/denominator )
  return( vip )
}
