
#' Bulk data from Simulation Study 3
#'
#' This observed real bulk Hi-C data correspond to the long arm of chromosome 16 of the K562 bulk A dataset with 200 kb resolution (GEO accession: GSM2109887).
#'
#' @format A symmetric count matrix with dimension 220 by 220. The proportion of observed zeros is 79.2%.
#'
#' @source https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2109887, https://github.com/osu-stat-gen/RWS-RWR
"sim3.bulk"


#' Knight-Ruiz (KR) normalized bulk data from Simulation Study 3
#'
#' This data is obtained by performing KR-normalization on the `sim3.bulk` data.
#'
#' @format A symmetric square matrix, with row sums (and column sums) all equal to 1.
#'
#' @source https://github.com/osu-stat-gen/RWS-RWR
"sim3.bulk.KR"


#' One subsampled data from Simulation Study 3
#'
#' This data is obtained by subsampling 10% of the observed counts in `sim3.bulk`.
#'
#' @format A symmetric count matrix with dimension 220 by 220. The proportion of observed zeros is 95.5%, which is about the same as a high-quality single-cell Hi-C matrix.
#'
#' @source https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2109887, https://github.com/osu-stat-gen/RWS-RWR
"sim3.sub"


#' Knight-Ruiz (KR) normalized subsampled data from Simulation Study 3
#'
#' This data is obtained by performing KR-normalization on the `sim3.sub` data.
#'
#' @format A symmetric square matrix, with row sums (and column sums) all equal to 1.
#'
#' @source https://github.com/osu-stat-gen/RWS-RWR
"sim3.sub.KR"

