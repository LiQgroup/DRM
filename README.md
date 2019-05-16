# DRM
Distance Regression Model

## Description
Distance-based regression model is a commonly used statistical analysis strategy to identify feature patterns in genomics, genetics and other application areas. It is specially useful for high-dimensional biological or genetic data with a proper distance or similarity measure available. This model brings in a pseudo F test statistic and its signiﬁcance is computed via permutation procedure, which is computationally demanding and might restrict its application in large-scale studies. We ﬁnd that when the sample size is large, the denominator of the pseudo F test statistic is almost invariant. Therefore, the numerator is a good approximation of the F test statistic. We standardize the numerator and derive its asymptotic distribution to obtain the p-value when the sample size is large. 

The package “DRM” contains two functions that can calculate p-values in the distance regression model, including drm.pf (pseudo F test statistic) and fast.pf (standard error difference, calculated using asymptotic distribution or permutation procedure).

## Dependencies on R packages
* MASS

## Installation
Install the R package “DRM” using the R package “devtools”.

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

2.	Overview of the package “DRM”

```
help(package=DRM)
```

3.	A toy example for the function “simreg”

The function “simreg” calculates the p-value for the pseudo F test statistic using permutation procedure. Here is a toy example.

```
# Generate the matrix which combines the covariate matrix and the matrix of interest.
> mu.x <- rep(2, 10)
> sigma.x <- diag(rep(1,10))
> n <- 30
> X <- mvrnorm(n, mu.x, sigma.x)
```
