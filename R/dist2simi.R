##' Transform a distance matrix to a similarity matrix.
##'
##' Available transformation methods are:
##' \tabular{llll}{
##' \code{"m"} \tab \tab \tab \cr
##'  \tab \tab \tab \eqn{s_{ij} = -0.5*d_{ij}^2} \cr
##' \code{"e"} \tab \tab \tab \cr   
##'  \tab \tab \tab \eqn{s_{ij} = exp(-0.5*d_{ij}^2)} \cr
##' \code{"es"} \tab \tab \tab \cr
##'  \tab \tab \tab   \eqn{s_{ij} = exp(-0.5*d_{ij}^2/mean(D^2))} \cr
##' \code{"r"} \tab \tab \tab \cr
##'  \tab \tab \tab    \eqn{s_{ij} = 1/(1+d_{ij})} \cr
##' \code{"rs"} \tab \tab \tab \cr
##'  \tab \tab \tab   \eqn{s_{ij} = 1/(1+d_{ij}/mean(D))} \cr
##' }
##' @title Transform a distance matrix to a similarity matrix
##' @param D A distance matrix.
##' @param method The method used to conduct transformation. This must
##' be one of \code{"m"}, \code{"e"}, \code{"es"}, \code{"r"},
##' \code{"rs"}. The default is \code{"m"}. See details.
##' @return A similarity matrix.
##' @author Jinjuan Wang.
##' @export
##' @examples
##' library(MASS)
##' mu.x = rep(2, 10)
##' sigma.x = diag(rep(1,10))
##' n = 30
##' X = mvrnorm(n, mu.x, sigma.x)
##' B = matrix(rnorm(10*100, mean = 3, sd = 3), nrow = 10)
##' Y = X %*% B
##' D = dist(Y)
##' S = dist2simi(D, "m")

dist2simi <- function(D, method = "m")
{
  if(method == "m"){
    S = -0.5*D^2
  }else if(method == "e"){
    S = exp(-0.5*D^2)
  }else if(method == "es"){
    S = exp(-0.5*D^2/mean(D^2))
  }else if(method == "r"){
    S = 1/(1+D)
  }else if(method == "rs"){
    S = 1/(1+D/mean(D))
  }

  return(S)

}
