
<!-- README.md is generated from README.Rmd. Please edit that file -->

# RW

<!-- badges: start -->
<!-- badges: end -->

The goal of RW package is to investigate whether random-walk based
methods — random walk with limited number of steps (RWS) and random walk
with restart (RWR) — can improve Hi-C data quality, especially for the
detection of Topologically Associating Domains (TADs). In this package,
we introduced three simulation studies, along with two real datasets, to
conduct the investigation at different sparsity levels and TAD clarity
levels.

## Installation

You can install the RW package from GitHub:

``` r
# If not already installed:
# install.packages("devtools")

# library(devtools)
# devtools::install_github("https://github.com/osu-stat-gen/RW")
library(RW)
```

## Example

### Data Preparation: KR normalization

In order to perform random-walk based methods (RWS and RWR), we need to
transform the Hi-C contact count matrix into a transition probability
matrix (TPM). Here we apply Knight-Ruiz (KR) algorithm to turn a contact
matrix (symmetric and all elements being non-negative) into a symmetric
TPM.

``` r
# Hi-C contact count for long arm of chromosome 22 of a bulk hESC dataset
data("hESC.bulk")
# Perform KR normalization
hESC.bulk.KR <- KRnorm2(hESC.bulk)
#> Cols/Rows being zeros: 33 of them
#>   41  42  43  44  45  82  83  84  85  86  87  88  89  90  99  110  111  112  113  114  115  116  117  119  120  129  199  743  840  854  855  856  857

print(dim(hESC.bulk.KR))  # same as the dimension of hESC.bulk, which is 857 by 857. 
#> [1] 857 857
```

### Perform RWS or RWR

Once we obtain the TPM based on Hi-C data, we can perform either random
walk with limited number of steps (RWS) or random walk with restart
(RWR). In our package, the two functions are called `RWS()` and `RWR()`
respectively.

`RWS()` takes the TPM $P$ (say, a KR-normalized contact matrix) as the
first argument, the number of steps $k$ as the second argument, and
computes the $k$-step transition matrix $P^{k}$.

``` r
# Perform RWS with 3 steps
hESC.bulk.RWS3s <- RWS(hESC.bulk.KR, step=3)
```

`RWR()` calculates the limiting-state RWR matrix, given the TPM $P$ and
restart probability $\alpha$. The computation is done through
fixed-point iteration, which stops when the Frobenius norm of the
difference between the last two iterations is less than a small
tolerance `tol` or when `max_iter` is reached.

``` r
# Perform RWR with alpha=0.1 
hESC.bulk.RWR0.1 <- RWR(hESC.bulk.KR, alpha=0.1, max_iter=300, tol=1e-06) 
#> Iterations: 49; last change (Frobenius norm): 8.481e-07
```

For more details, please see the vignette.
