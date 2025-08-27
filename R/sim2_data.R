
#' Simulated \eqn{E} data from Simulation Study 2
#'
#' \eqn{E} is a simulated Hi-C data with `N=1000` bins and `k=20` TADs. The domain endpoints are saved as a vector called `sim2.true.bound.E`.
#'
#' @format A symmetric square matrix with 20 diagonal blocks and many random background noise.
#'
#' @source https://github.com/osu-stat-gen/RWS-RWR
"sim2.E"


#' Knight-Ruiz (KR) normalized \eqn{E} data from Simulation Study 2
#'
#' This data is obtained by performing KR-normalization on the `sim2.E` data.
#'
#' @format A symmetric square matrix, with row sums (and column sums) all equal to 1.
#'
#' @source https://github.com/osu-stat-gen/RWS-RWR
"sim2.E.KR"


#' Domain boundary points of the \eqn{E} data from Simulation Study 2
#'
#' The `N=1000` bins in the `sim2.E` data are divided consecutively into `k=20` TADs. The 20 endpoints are listed here.
#'
#' @format An integer vector marking the endpoints of the 20 TADs.
#'
#' @source https://github.com/osu-stat-gen/RWS-RWR
"sim2.true.bound.E"

