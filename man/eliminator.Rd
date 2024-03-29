% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eliminator.R
\name{eliminator}
\alias{eliminator}
\title{The Eliminator algorithm}
\usage{
eliminator(
  y,
  X,
  reg = 0.5,
  prior = NULL,
  max.dim = NULL,
  frac = 0.25,
  vip.lim = 1,
  n.seg = 10,
  verbose = TRUE
)
}
\arguments{
\item{y}{Vector of responses, a factor of exact 2 levels.}

\item{X}{Matrix of predictor values.}

\item{reg}{The regularization parameter, see below.}

\item{prior}{Vector of prior probabilities, one value for each factor level in \code{y}.}

\item{max.dim}{Integer, the maximum number of dimensions to consider.}

\item{frac}{Fraction of unimportant variables to eliminate in each iteration (default is 0.25).}

\item{vip.lim}{The threshold for the VIP criterion (default is 1.0).}

\item{n.seg}{Integer, the number of cross-validation segments (default is 10).}

\item{verbose}{Logical, turns on/off output during computations.}
}
\value{
A \code{list} with two matrices, \code{Elimination} and \code{Selected}.

\code{Elimination} has one row for each iteration, with accuracy results after each. The first column is
the number of variables left, the second is the fraction of correctly classified samples, and the third
is the p-value of the McNemar-test, see above.
You will typically look down this matrix for iterations where you have the maximum accuracy. Then,
continue down the matrix from that row, until there is a significant drop in accuracy (column 2) and a
corresponding small p-value (column 3), and choose the result just
before this as your final selection.

In \code{Selected} you also find one row for each iteration. This logical indicates which variables are
selected after each iteration.
Thus, if cell \code{[i,j]} is \code{TRUE} variable \code{j} is still included after iteration \code{i}.
}
\description{
Variable elimination using pda.
}
\details{
This is a variable selection (elimination) algorithm based on the \code{\link{pda}} method
where the response (\code{y}) is a factor with 2 levels, i.e. a two-class problem. The restriction
to a two-class problem comes from the use of the VIP criterion. The idea is a slight modification
of the one described in Mehmood et al, (2011).

The algorithm starts out by using all variables, fitting a \code{\link{pda}} model, using the
\code{\link{pdaDim}} function to estimate the optimal dimensionality for the PLS step after each
elimination. The regularization parameter \code{reg} and the number of
cross-validation segments \code{n.seg} are needed for this, see \code{\link{pdaDim}} for details.

Based on the fitted model, the VIP-criterion is computed for all variables, see \code{\link{vipCriterion}}
for details. Variables are ranked by this criterion, and those with a VIP-score below \code{vip.lim}
are considered unimportant. The fraction \code{frac} of the unimportant variables are eliminated.
If no unimportant variables are left, the algorithm terminates. It will mostly continue until
only 1 variable is left, and then terminate.

After each elimination, the classification accuracy, i.e. the fraction of correctly classified samples, is
computed based on the cross-validation of the \code{\link{pdaDim}} results. The reduced model accuracy is
tested against the maximum so far, using the \code{\link{mcnemar.test}}. The p-value of this test is
reported for each elimination step, indicating if
the reduced model has a significantly poorer accuracy than the maximum so far. The idea is to keep
on eliminating variables until the accuracy becomes significantly poorer than the maximum.

The argument \code{max.dim} allows you to specify the maximum number of PLS components to consider. The
optimal dimension found by \code{\link{pdaDim}} is between 1 and \code{max.dim}. If it turns out identical
to your \code{max.dim}, increase its value slightly and re-run.
}
\examples{
data(microbiome)
y <- microbiome[1:40, 1]
X <- as.matrix(microbiome[1:40, -1])
lst <- eliminator(y, X, max.dim = 10)
print(lst$Elimination)
# Seems like iteration 23 is the place to stop, since a significant drop
# in performance is found at iteration 24.
# There are 2 variables left in iteration 23. These are variables
print(which(lst$Selected[23,]))

}
\references{
Mehmood, T, Martens, H, Saeboe, S, Warringer, J, Snipen, L (2011). A Partial Least Squares based algorithm
for parsimonious variable selection. Algorithms for Molecular Biology, 6:27.
}
\seealso{
\code{\link{vipCriterion}}, \code{\link{pdaDim}}.
}
\author{
Lars Snipen.
}
