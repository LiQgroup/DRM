##' Calculate the \emph{p}-value for the pseudo \emph{F} test
##' statistic using the permutation procedures.
##'
##' \code{x.mat[, null.space]} indicates the covariate matrix and 
##' \code{x.mat[, -null.space]} indicates the predictors of interest.
##' @title Pseudo F test in distance-based regression
##' @param simi.mat A similarity matrix among the subjects.
##' @param null.space A numeric vector containing the indices of those
##' columns in \code{x.mat} corresponding to the predictors which are
##' to be adjusted.
##' @param x.mat The matrix which combines the covariate matrix and
##' the matrix of interest.
##' @param permute A logical value indicating whether the
##' \emph{p}-value is calculated based on the permutation procedures. The
##' default is \code{TRUE}.
##' @param n.monterCarlo The repeat number of the permutation
##' procedures. The default is \code{1000}.
##' @param seed The seed of the random number generator for the
##' permutation procedures. The default is \code{NULL}.
##' @return The \emph{p}-value for the pseudo \emph{F} test statistic.
##' @author Qizhai Li.
##' @export
##' @examples
##' library(MASS)
##' mu.x <- rep(2, 10)
##' sigma.x <- diag(rep(1,10))
##' B <- matrix(rnorm(10*100, mean = 3, sd = 3), nrow = 10)
##' n <- 30
##' X <- mvrnorm(n, mu.x, sigma.x)
##' Y <- X %*% B
##' S <- Y %*% t(Y)
##' drm.pf(S, 1:5, X)
'drm.pf' <- function(simi.mat, null.space, x.mat, permute=TRUE, n.monterCarlo=1000, seed=NULL)
{
  if (!is.null(seed))
  {
    set.seed(seed)
  }

  x1 <- x.mat[,null.space]

  if (length(null.space)==1)
  {
    x1 <- matrix(x1, ncol=1)
  }

  x.hat <- x.mat %*% solve(t(x.mat)%*%x.mat) %*% t(x.mat)

  x1.hat <- x1 %*% solve(t(x1)%*%x1) %*% t(x1)

  n <- nrow(simi.mat)

  I.n <- diag(n)

  cent <- I.n - matrix(1,nrow=n,ncol=n)/n

  i.x  <- I.n - x.hat

  i.x1 <- I.n - x1.hat

  Q <- i.x1 %*% cent %*% simi.mat %*% cent %*% i.x1

  alter.hat <- x.hat - x1.hat

  F.obs <- sum(alter.hat*Q) / sum(i.x*Q)

  U <- 1:n

  F.star <- rep(NA, n.monterCarlo)

  for (i in 1:n.monterCarlo)
  {
    id.sam <- sample(U, replace=!permute)

    Q.star <- Q[id.sam, id.sam]

    F.star[i] <- sum(alter.hat*Q.star) / sum(i.x*Q.star)
  }

  sum(F.star >= F.obs)/n.monterCarlo
}
