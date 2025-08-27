
#' Choose the process region for TopDom (internal)
#' @noRd
#' @keywords internal
Which.process.region_star <- function (rmv.idx, n_bins, min.size = 3L){
  gap.idx <- rmv.idx
  proc.regions <- data.frame(start = numeric(0), end = numeric(0),
                             stringsAsFactors = FALSE)
  proc.set <- setdiff(seq_len(n_bins), gap.idx)
  n_proc.set <- length(proc.set)
  i <- 1
  while (i <= n_proc.set) {
    start <- proc.set[i]
    j <- i + 1
    while (j <= n_proc.set) {
      if (proc.set[j] - proc.set[j - 1] <= 1) {
        j <- j + 1
      }
      else {
        proc.regions <- rbind(proc.regions, c(start = start,
                                              end = proc.set[j - 1]))
        i <- j
        break
      }
    }

    if (j > n_proc.set) {
      proc.regions <- rbind(proc.regions, c(start = start,
                                            end = proc.set[j - 1]))
      break
    }
  }
  colnames(proc.regions) <- c("start", "end")
  proc.regions <- proc.regions[abs(proc.regions[, "end"] -
                                     proc.regions[, "start"]) >= min.size, ]
  proc.regions
}


#' Convert bins to domains for TopDom (internal)
#' @noRd
#' @keywords internal
Convert.Bin.To.Domain.TMP_star <- function (bins, signal.idx, gap.idx, pvalues = NULL, pvalue.cut = NULL){
  n_bins <- nrow(bins)
  ret <- data.frame(chr = character(0), from.id = numeric(0),
                    from.coord = numeric(0), to.id = numeric(0), to.coord = numeric(0),
                    tag = character(0), size = numeric(0), stringsAsFactors = FALSE)
  levels(x = ret[, "tag"]) <- c("domain", "gap",
                                "boundary")
  rmv.idx <- setdiff(seq_len(n_bins), gap.idx)
  proc.region <- Which.process.region_star(rmv.idx, n_bins = n_bins,
                                           min.size = 0L)
  from.coord <- bins[proc.region[, "start"], "from.coord"]
  n_procs <- nrow(proc.region)
  zeros <- double(length = n_procs)
  gap <- data.frame(chr = rep(bins[1, "chr"], times = n_procs),
                    from.id = zeros, from.coord = from.coord, to.id = zeros,
                    to.coord = zeros, tag = rep("gap", times = n_procs),
                    size = zeros, stringsAsFactors = FALSE)
  rmv.idx <- union(signal.idx, gap.idx)
  proc.region <- Which.process.region_star(rmv.idx, n_bins = n_bins,
                                           min.size = 0L)
  n_procs <- nrow(proc.region)
  from.coord <- bins[proc.region[, "start"], "from.coord"]
  zeros <- double(length = n_procs)
  domain <- data.frame(chr = rep(bins[1, "chr"], times = n_procs),
                       from.id = zeros, from.coord = from.coord, to.id = zeros,
                       to.coord = zeros, tag = rep("domain", times = n_procs),
                       size = zeros, stringsAsFactors = FALSE)
  rmv.idx <- setdiff(seq_len(n_bins), signal.idx)
  proc.region <- as.data.frame(Which.process.region_star(rmv.idx,
                                                         n_bins = n_bins, min.size = 1L))
  n_procs <- nrow(proc.region)
  if (n_procs > 0) {
    from.coord <- bins[proc.region[, "start"] + 1,
                       "from.coord"]
    zeros <- double(length = n_procs)
    boundary <- data.frame(chr = rep(bins[1, "chr"],
                                     times = n_procs), from.id = zeros, from.coord = from.coord,
                           to.id = zeros, to.coord = zeros, tag = rep("boundary",
                                                                      times = n_procs), size = zeros, stringsAsFactors = FALSE)
    ret <- rbind(ret, boundary)
  }
  if (nrow(domain) == 0L) {
    ret <- gap
  }
  else {
    ret <- rbind(gap, domain)
    ret <- ret[order(ret[, 3]), ]
    ret[, "to.coord"] <- c(ret[2:nrow(ret), "from.coord"],
                           bins[n_bins, "to.coord"])
    ret[, "from.id"] <- match(ret[, "from.coord"],
                              table = bins[, "from.coord"])
    ret[, "to.id"] <- match(ret[, "to.coord"],
                            table = bins[, "to.coord"])
    ret[, "size"] <- ret[, "to.coord"] - ret[,
                                             "from.coord"]
    if (!is.null(pvalues) && !is.null(pvalue.cut)) {
      for (i in seq_len(nrow(ret))) {
        if (ret[i, "tag"] == "domain") {
          domain.bins.idx <- ret[i, "from.id"]:ret[i,
                                                   "to.id"]
          p.value.constr <- which(pvalues[domain.bins.idx] <
                                    pvalue.cut)
          if (length(domain.bins.idx) == length(p.value.constr))
            ret[i, "tag"] <- "boundary"
        }
      }
    }
  }
  ret
}


#' Get p-value of each bin for TopDom (internal)
#' @noRd
#' @keywords internal
#' @importFrom rstatix wilcox_test
Get.Pvalue_star <- function (matrix.data, size, scale = 1){
  n_bins <- nrow(matrix.data)
  pvalue <- rep(1, times = n_bins)
  for (i in seq_len(n_bins - 1)) {
    dia <- as.vector(TopDom:::Get.Diamond.Matrix2(matrix.data, i = i, size = size))
    ups <- as.vector(TopDom:::Get.Upstream.Triangle(matrix.data, i = i, size = size))
    downs <- as.vector(TopDom:::Get.Downstream.Triangle(matrix.data, i = i, size = size))

    test.data <- data.frame(y=c(dia, ups, downs), x=factor(rep(c("dia","tri"), c(length(dia), length(ups)+length(downs)))))
    coin.wil.test <- wilcox_test(y ~ x, data = test.data, exact=TRUE, alternative = "less")  # conf.int = FALSE
    pvalue[i] <- pvalue(coin.wil.test)
  }
  pvalue[is.na(pvalue)] <- 1
  pvalue
}


#' @keywords internal
#' @noRd
#' @importFrom utils getFromNamespace
td_Get.Diamond.Matrix <- function(...) {
  getFromNamespace("Get.Diamond.Matrix", "TopDom")(...)
}

td_Which.Gap.Region2 <- function(...) {
  getFromNamespace("Which.Gap.Region2", "TopDom")(...)
}

td_Detect.Local.Extreme <- function(...) {
  getFromNamespace("Detect.Local.Extreme", "TopDom")(...)
}





# Extract the upper triangle of the contact matrix by row (helper function for Subsampling, internal)
#' @noRd
#' @keywords internal
matrixtovec <- function(M){
  vec <- c()
  m <- dim(M)[1]
  for(i in 1:m){
    vec <- c(vec, M[i, i:m])
  }
  return(vec)
}

# Transform a vector of length n*(n+1)/2 back to a symmetric m*m matrix (helper function for Subsampling, internal)
#' @noRd
#' @keywords internal
vectomatrix <- function(vec){
  m <- floor(sqrt(2*length(vec)))
  M <- matrix(0, nrow=m, ncol=m)
  index <- 0
  for(i in 1:m){
    M[i,i:m] <- vec[(index+1):(index+m-i+1)]
    index <- (index+m-i+1)
  }
  diag(M) <- 0.5*diag(M)
  M <- (M+t(M))
  return(M)
}



# Function that takes the heatmap and boundary points as input, and returns the heatmap with lines marking the blocks #
#' @noRd
#' @keywords internal
add_bound <- function(normplot, bound, col="red"){
  plot1 <- normplot
  bound <- c(0, bound[bound>0])
  nb <- length(bound)
  if(nb<=2){
    segment = data.frame(
      x = c(rep(0.5,3), bound[2]+0.5),
      y = c(bound[2]+0.5, 0.5, rep(bound[2]+0.5, 2)),
      xend = c(rep(bound[2]+0.5, 2), 0.5, bound[2]+0.5),
      yend = c(bound[2]+0.5, rep(0.5,3))
    )
    plot1 <- plot1 + annotate("segment", x = segment$x, xend = segment$xend, y = segment$y, yend = segment$yend, color=col, linewidth=0.1)
    return(plot1)
  } else{
    segment = data.frame(
      x = c(bound[2:(nb-1)]+0.5, bound[1:(nb-2)]+0.5),
      y = c(bound[nb]+0.5-bound[1:(nb-2)], bound[nb]+0.5-bound[2:(nb-1)]),
      xend = c(bound[2:(nb-1)]+0.5, bound[3:nb]+0.5),
      yend = c(bound[nb]+0.5-bound[3:nb], bound[nb]+0.5-bound[2:(nb-1)])
    )
    plot1 <- plot1 + annotate("segment", x = segment$x, xend = segment$xend, y = segment$y, yend = segment$yend, color=col, linewidth=0.1)
    return(plot1)
  }
}


