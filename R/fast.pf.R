##' Calculate the \emph{p}-value for the standard error difference
##' test using the asymptotic distribution or the permutation
##' procedures.
##'
##' If \code{permutation} is \code{FALSE}, the \emph{p}-value is
##' calculated using the standard normal distribution.

##' @title Standard error difference in distance-based regression
##' @param X1 A matrix containing the predictors to be adjusted.
##' @param X2 A matrix containing the predictors of interest.
##' @param S0 A similarity matrix.
##' @param permutation A logical value indicating whether the
##' permutation procedure is used. If \code{FALSE}, the \emph{p}-value
##' is calculated by the asymptotic distribution. The default is
##' \code{FALSE}.
##' @param perNum The repeat number of the permutation. The default is
##' \code{1000}.
##' @return The \emph{p}-value for the standard error difference test.
##' @author Jinjuan Wang.
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
##' X1 <- X[, 1:5]
##' X2 <- X[, 6:10]
##' fast.pf(X1, X2, S, permutation = FALSE)
'fast.pf' <- function(X1, X2, S0, permutation = FALSE, perNum = 1000)
{
  X = cbind(X1, X2)

  if(is.matrix(X1) == TRUE){
    n = dim(X1)[1]
    k1 = dim(X1)[2]
  }else{
    n = length(X1)
    k1 = 1
  }

  if(is.matrix(X2) == TRUE){
    k2 = dim(X2)[2]
  }else{
    k2 = 1
  }

  k = k1 + k2

  H = diag(rep(1,n)) - matrix(1/n, nrow = n, ncol = n)
  S = H%*%S0%*%H

  if(permutation == FALSE){

    eigS = eigen(S)
    eigVal = eigS$values

    minVal = min(eigVal)

    if(minVal < -0.001){

      print("The similarity matrix is not semi-definite.")

      pvalue1 = NA

    }else{
      u = sum(eigVal>0.00005)
      eigMat = eigS$vectors[,1:u]
      newY = eigMat%*%diag(sqrt(eigVal[1:u]))
      # newS = newY%*%t(newY)

      Hx = X%*%ginv(crossprod(X))%*%t(X)
      Hx1 = X1%*%ginv(crossprod(X1))%*%t(X1)

      diffH1 = Hx - Hx1
      diffH2 = diag(rep(1,n)) - Hx

      del = (t(newY)%*%diffH2%*%newY)/(n-k)

      v1 = sum(diffH1*S)
      mu1 = sum(diag(diffH1))*sum(diag(del))
      sigma1 = 2*sum(diffH1^2)*((sum(del^2)-((sum(diag(del)))^2)/(n-k))*(n-k)^2/((n-k+2)*(n-k-1)))
      w1 = (v1-mu1)/sqrt(sigma1)

      pvalue1 = 1 - pnorm(w1)
    }



  }else{

    Hx = X%*%ginv(crossprod(X))%*%t(X)
    Hx1 = X1%*%ginv(crossprod(X1))%*%t(X1)

    I.n = diag(n)
    i.x1 = I.n - Hx1
    Snew1 = i.x1 %*% S %*% i.x1

    diffH1 = Hx - Hx1

    U = 1:n

    tVec = rep(NA, perNum)
    for (i in 1:perNum){
      id.sam = sample(U, replace= FALSE)

      Snew.star = Snew1[id.sam, id.sam]

      tVec[i] = sum(diffH1*Snew.star)
    }

    t1 = sum(diffH1*Snew1)
    pvalue1 = sum(tVec>t1)/perNum
  }

  return(pvalue1)

}

