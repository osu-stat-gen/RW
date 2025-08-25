
#' Corrected function of TopDom for TAD detection
#'
#' `TopDom_star()` computes bin-level contact signals, detects local extrema as
#' candidate TAD boundaries, optionally performs statistical filtering, and
#' converts boundaries to domains (i.e., TADs).
#'
#' @importFrom utils read.table write.table
#' @description
#' The function follows three steps:
#' \enumerate{
#'   \item Compute per-bin contact frequency using a diamond window of \code{window.size}.
#'   \item Detect gap regions and local extrema (candidate boundaries).
#'   \item (Optional) Statistical filtering via diagonal-wise scaling and
#'         Wilcoxon rank-sum p-values to reduce false positives.
#' }
#'
#' @param matrix.file Either (a) an object of class \code{"TopDomData"} with
#'   components \code{bins} (data.frame) and \code{counts} (numeric matrix), or
#'   (b) a character path to a TopDom-format text file whose first 3–4 columns
#'   describe bins \code{(chr, from.coord, to.coord[, id])} followed by an
#'   \eqn{n \times n} symmetric contact matrix, where `n` is the number of bins.
#' @param window.size Integer (>=1). Diamond half-size used for local mean
#'   calculations and subsequent statistics.
#' @param outFile Optional character prefix. If supplied, three files are
#'   written: \code{<prefix>.binSignal}, \code{<prefix>.domain}, and
#'   \code{<prefix>.bed}.
#' @param statFilter Logical; if \code{TRUE} (default), perform Step 3 statistical
#'   filtering of candidate boundaries.
#'
#' @return A list with three components:
#' \describe{
#'   \item{binSignal}{\code{data.frame} with per-bin signals and p-values
#'         (\code{chr}, \code{from.coord}, \code{to.coord}, \code{local.ext},
#'         \code{mean.cf}, \code{pvalue}).}
#'   \item{domain}{\code{data.frame} of called domains.}
#'   \item{bed}{\code{data.frame} with BED fields
#'         \code{chrom}, \code{chromStart}, \code{chromEnd}, \code{name}.}
#' }
#'
#' @details
#' This function is a modification the `TopDom()` function from the cran package built
#' by the original authors. The process regions are adjusted so that TopDom can handle
#' the existence of bins that never interact with other bins without excluding them
#' completely from the analysis.
#'
#' @references
#' Shin, H., Shi, Y., Dai, C., Tjong, H., Gong, K., Alber, F., & Zhou, X. J. (2016).
#' TopDom: an efficient and deterministic method for identifying topological domains in genomes.
#' \emph{Nucleic Acids Research}, 44(7), e70. \doi{10.1093/nar/gkv1505}
#' PMCID: \href{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4838359}{PMC4838359};
#' PMID: \href{https://pubmed.ncbi.nlm.nih.gov/26704975/}{26704975}.
#'
#' @seealso \code{\link[TopDom]{TopDom}}
#'
#' @examples
#' \dontrun{
#' data(sim1.P.obs)
#' N <- nrow(sim1.P.obs)
#' P.obs.input <- data.frame(rep(paste0("chr", as.character(1)), N),
#'                           (0:(N-1))*40000,
#'                           (1:N)*40000,
#'                           sim1.P.obs)
#' data.table::fwrite(P.obs.input, file="P_obs_TopDom.txt", col.names=FALSE, sep=" ")
#' P_obs_TopDom_data <- TopDom::readHiC("P_obs_TopDom.txt", binSize = 40000)
#' tryCatch({
#'   P_obs_TopDom <- TopDom_star(P_obs_TopDom_data, window.size=5)
#'   topdom.bound.P.obs <- P_obs_TopDom$domain[,"to.id"]},
#'   error = function(c) "error"
#' )
#' if(exists("P_obs_TopDom") == FALSE){
#'   topdom.bound.P.obs <- c(N)
#' }
#' }
#' @export
TopDom_star <- function(matrix.file, window.size, outFile=NULL, statFilter=TRUE){
  if (inherits(matrix.file, "TopDomData")) {
    bins <- matrix.file$bins
    matrix.data <- matrix.file$counts
    n_bins <- as.double(nrow(bins))
    mean.cf <- rep(0, times = n_bins)
    pvalue <- rep(1.0, times = n_bins)
    local.ext <- rep(-0.5, times = n_bins)
  } else {
    print("#########################################################################")
    print("Step 0 : File Read ")
    print("#########################################################################")
    window.size = as.numeric(window.size)
    matdf <- read.table(matrix.file, header=F)

    if( ncol(matdf) - nrow(matdf) == 3) {
      colnames(matdf) <- c("chr", "from.coord", "to.coord")
    } else if( ncol(matdf) - nrow(matdf) ==4 ) {
      colnames(matdf) <- c("id", "chr", "from.coord", "to.coord")
    } else {
      print("Unknwon Type of matrix file")
      return(0)
    }
    n_bins = nrow(matdf)
    mean.cf <- rep(0, n_bins)
    pvalue <- rep(1, n_bins)

    local.ext = rep(-0.5, n_bins)

    bins <- data.frame(id=1:n_bins,
                       chr=matdf[, "chr"],
                       from.coord=matdf[, "from.coord"],
                       to.coord=matdf[, "to.coord"] )
    matrix.data <- as.matrix( matdf[, (ncol(matdf) - nrow(matdf)+1 ):ncol(matdf)] )

    print("-- Done!")
    print("Step 0 : Done !!")
  }


  print("#########################################################################")
  print("Step 1 : Generating binSignals by computing bin-level contact frequencies")
  print("#########################################################################")
  ptm <- proc.time()
  for(i in 1:n_bins)
  {
    diamond = td_Get.Diamond.Matrix(mat.data=matrix.data, i=i, size=window.size)
    mean.cf[i] = mean(diamond)
  }

  eltm = proc.time() - ptm
  print(paste("Step 1 Running Time : ", eltm[3]))
  print("Step 1 : Done !!")

  print("#########################################################################")
  print("Step 2 : Detect TD boundaries based on binSignals")
  print("#########################################################################")

  ptm = proc.time()
  #gap.idx = Which.Gap.Region(matrix.data=matrix.data)
  #gap.idx = Which.Gap.Region2(mean.cf)
  gap.idx = td_Which.Gap.Region2(matrix.data=matrix.data, window.size)

  proc.regions = Which.process.region_star(rmv.idx=gap.idx, n_bins=n_bins, min.size=3)

  #print(proc.regions)

  for( i in 1:nrow(proc.regions))
  {
    start = proc.regions[i, "start"]
    end = proc.regions[i, "end"]

    print(paste("Process Regions from ", start, "to", end))

    local.ext[start:end] = td_Detect.Local.Extreme(x=mean.cf[start:end])
  }

  eltm = proc.time() - ptm
  print(paste("Step 2 Running Time : ", eltm[3]))
  print("Step 2 : Done !!")

  if(statFilter)
  {
    print("#########################################################################")
    print("Step 3 : Statistical Filtering of false positive TD boundaries")
    print("#########################################################################")

    ptm = proc.time()
    print("-- Matrix Scaling....")
    n_bins <- as.double(n_bins)
    print(as.double(n_bins*n_bins))
    scale.matrix.data = matrix.data
    for( i in 1:(2*window.size) )
    {
      #diag(scale.matrix.data[, i:n_bins] ) = scale( diag( matrix.data[, i:n_bins] ) )
      scale.matrix.data[ seq(1+(n_bins*i), as.double(n_bins*n_bins), 1+n_bins) ] = scale( matrix.data[ seq(1+(n_bins*i), as.double(n_bins*n_bins), 1+n_bins) ] )
    }

    print("-- Compute p-values by Wilcox Ranksum Test")
    for( i in 1:nrow(proc.regions))
    {
      start = proc.regions[i, "start"]
      end = proc.regions[i, "end"]

      print(paste("Process Regions from ", start, "to", end))

      pvalue[start:end] <- Get.Pvalue_star(matrix.data=scale.matrix.data[start:end, start:end], size=window.size, scale=1)
    }
    print("-- Done!")

    print("-- Filtering False Positives")
    local.ext[intersect( union(which( local.ext==-1), which(local.ext==-1)), which(pvalue<0.05))] = -2
    local.ext[which(local.ext==-1)] = 0
    local.ext[which(local.ext==-2)] = -1
    print("-- Done!")

    eltm = proc.time() - ptm
    print(paste("Step 3 Running Time : ", eltm[3]))
    print("Step 3 : Done!")
  } else pvalue = 0

  domains = Convert.Bin.To.Domain.TMP_star(bins=bins,
                                           signal.idx=which(local.ext==-1),
                                           gap.idx=which(local.ext==-0.5),
                                           pvalues=pvalue,
                                           pvalue.cut=0.05)

  bins = cbind(bins,
               local.ext = local.ext,
               mean.cf = mean.cf,
               pvalue = pvalue)

  bedform = domains[, c("chr", "from.coord", "to.coord", "tag")]
  colnames(bedform) = c("chrom", "chromStart", "chromEnd", "name")

  if( !is.null(outFile) ) {
    print("#########################################################################")
    print("Writing Files")
    print("#########################################################################")

    outBinSignal =  paste(outFile, ".binSignal", sep="")
    print(paste("binSignal File :", outBinSignal) )
    write.table(bins, file=outBinSignal, quote=F, row.names=F, col.names=T, sep="\t")


    outDomain = paste(outFile, ".domain", sep="")
    print(paste("Domain File :", outDomain) )
    write.table( domains, file=outDomain, quote=F, row.names=F, col.names=T, sep="\t")

    outBed = paste(outFile, ".bed", sep="")
    print(paste("Bed File : ", outBed))
    write.table( bedform, file=outBed, quote=F, row.names=F, col.names=F, sep="\t")
  }

  print("Done!!")

  print("Job Complete !")
  return(list(binSignal=bins, domain=domains, bed=bedform))
}

