
#' Human Embryonic Stem Cells (hESC) bulk Hi-C data
#'
#' This observed real bulk Hi-C data correspond to the long arm of chromosome 22 of the hESC dataset with 40 kb resolution (GEO accession: GSE35156).
#'
#' @format A symmetric count matrix with dimension 857 by 857. The proportion of observed zeros is 73.8%.
#'
#' @source https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE35156, https://github.com/osu-stat-gen/RWS-RWR
"hESC.bulk"


#' Knight-Ruiz (KR) normalized human embryonic stem cells bulk Hi-C data
#'
#' This data is obtained by performing KR-normalization on the `hESC.bulk` data.
#'
#' @format A symmetric square matrix, with row sums (and column sums) all equal to 1.
#'
#' @source https://github.com/osu-stat-gen/RWS-RWR
"hESC.bulk.KR"

