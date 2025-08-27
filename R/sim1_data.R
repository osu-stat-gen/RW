
#' Simulated \eqn{P} data from Simulation Study 1
#'
#' \eqn{P} is a Hi-C data with `N=200` bins and `k=5` TADs. The domain sizes are specified as 50, 30, 20, 90, and 10, respectively.
#'
#' @format A symmetric square matrix with a clear 5-diagonal-block patterns and random background noise.
#'
#' @source https://github.com/osu-stat-gen/RWS-RWR
"sim1.P.obs"


#' Knight-Ruiz (KR) normalized \eqn{P} data from Simulation Study 1
#'
#' This data is obtained by performing KR-normalization on the `sim1.P.obs` data.
#'
#' @format A symmetric square matrix, with row sums (and column sums) all equal to 1.
#'
#' @source https://github.com/osu-stat-gen/RWS-RWR
"sim1.P.KR"

