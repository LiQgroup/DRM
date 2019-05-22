# DRM
Distance Regression Model

## Description
Distance-based regression model is a commonly used statistical analysis strategy to identify feature patterns in genomics, genetics and other application areas. It is specially useful for high-dimensional biological or genetic data with a proper distance or similarity measure available. This model brings in a pseudo F test statistic and its signiﬁcance is computed via permutation procedure, which is computationally demanding and might restrict its application in large-scale studies. We ﬁnd that when the sample size is large, the denominator of the pseudo F test statistic is almost invariant. Therefore, the numerator is a good approximation of the F test statistic. We standardize the numerator and derive its asymptotic distribution to obtain the p-value when the sample size is large. 

The package DRM contains two functions that can calculate p-values in the distance regression model, including drm.pf (pseudo F test statistic) and fast.pf (standard error difference, calculated using asymptotic distribution or permutation procedure). This package also contains a function which can transform a distance matrix to a similarity matrix.

## Dependencies on R packages
 * [MASS](https://cran.r-project.org/web/packages/MASS/index.html)

## Installation
Install the R package DRM using the R package [devtools](https://github.com/r-lib/devtools).

```
> install.packages("devtools")
> library(devtools) 
> install_github("LiQgroup/DRM")
```

## Short tutorial
1.	Import package

```
library(DRM)
```

2.	Overview of the package DRM

```
help(package=DRM)
```

3.	A toy example for the function “drm.pf”

The function “drm.pf” calculates the p-value for the pseudo F test statistic using permutation procedure. Here is a toy example.

```
# Generate the matrix which combines the covariate matrix and the matrix of interest.
> mu.x <- rep(2, 10)
> sigma.x <- diag(rep(1,10))
> n <- 30
> X <- mvrnorm(n, mu.x, sigma.x)

# Generate the similarity matrix.
> B <- matrix(rnorm(10*100, mean = 3, sd = 3), nrow = 10)
> Y <- X %*% B
> S <- Y %*% t(Y)

# Calculate the p-value for the pseudo F test statistic using permutation procedure.
> drm.pf(S, 1:5, X)
```
4.	A toy example for the function “fast.pf”

The function “fast.pf” calculates the p-value for the standard error difference test using the asymptotic distribution or the permutation procedure. Here is a toy example.

```
# Generate the matrix X1 which contains the predictors to be adjusted and the matrix X2 which contains the predictors of interest.
> mu.x <- rep(2, 10)
> sigma.x <- diag(rep(1,10))
> n <- 30
> X <- mvrnorm(n, mu.x, sigma.x)
> X1 <- X[, 1:5]
> X2 <- X[, 6:10]

# Generate the similarity matrix.
> B <- matrix(rnorm(10*100, mean = 3, sd = 3), nrow = 10)
> Y <- X %*% B
> S <- Y %*% t(Y)

# Calculate the p-value for the standard error difference test using the asymptotic distribution.
> fast.pf(X1, X2, S, permutation = FALSE)

# Calculate the p-value for the standard error difference test using the permutation procedure.
> fast.pf(X1, X2, S, permutation = TRUE)
```
5.	A toy example for the function “dist2simi”

The function “dist2simi” transform a distance matrix to a similarity matrix. Here is a toy example.

```
library(MASS)
# Generate a distance matrix.
mu.x = rep(2, 10)
sigma.x = diag(rep(1,10))
n = 30
X = mvrnorm(n, mu.x, sigma.x)
B = matrix(rnorm(10*100, mean = 3, sd = 3), nrow = 10)
Y = X %*% B
D = dist(Y)
# Transform the distance matrix to a similarity matrix.
S = dist2simi(D, "m")
```
## Citation
Jinjuan Wang, Wei Zhang, Wenjun Xiong, Lin Wang, and Qizhai Li. A quick and efﬁcient p-value calculation for distance-based regression model. Bioinformatics (submitted).
