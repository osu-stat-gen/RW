## Corrected TopDom ##
library(readr)
library(data.table)
library(scales)
library(TopDom)
library(coin)
#source("C:/R/R-4.1.0/library/TopDom/v0.03/fix_TopDom.R")




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




Get.Pvalue_star <- function (matrix.data, size, scale = 1){
  n_bins <- nrow(matrix.data)
  pvalue <- rep(1, times = n_bins)
  for (i in seq_len(n_bins - 1)) {
    dia <- as.vector(TopDom:::Get.Diamond.Matrix2(matrix.data, i = i, size = size))
    ups <- as.vector(TopDom:::Get.Upstream.Triangle(matrix.data, i = i, size = size))
    downs <- as.vector(TopDom:::Get.Downstream.Triangle(matrix.data, i = i, size = size))
    
    test.data <- data.frame(y=c(dia, ups, downs), x=factor(rep(c("dia","tri"), c(length(dia), length(ups)+length(downs)))))
    coin.wil.test <- wilcox_test(y ~ x, data = test.data, distribution = "exact", alternative = "less", conf.int = FALSE)
    pvalue[i] <- pvalue(coin.wil.test)
  }
  pvalue[is.na(pvalue)] <- 1
  pvalue
  
}




TopDom_star <- function( matrix.file, window.size, outFile=NULL, statFilter=T){
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
    diamond = TopDom:::Get.Diamond.Matrix(mat.data=matrix.data, i=i, size=window.size)
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
  gap.idx = TopDom:::Which.Gap.Region2(matrix.data=matrix.data, window.size)
  
  proc.regions = Which.process.region_star(rmv.idx=gap.idx, n_bins=n_bins, min.size=3)
  
  #print(proc.regions)
  
  for( i in 1:nrow(proc.regions))
  {
    start = proc.regions[i, "start"]
    end = proc.regions[i, "end"]
    
    print(paste("Process Regions from ", start, "to", end))
    
    local.ext[start:end] = TopDom:::Detect.Local.Extreme(x=mean.cf[start:end])
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


# 
# # An example 
# 
# data <- readHiC("K562_Bulk-B_40kb_chr21_TopDom.txt", binSize = 40000)
# tds_5 <- TopDom_0.0.3(data, window.size=5)
# tds_5
# tds_5_star <- TopDom_star(data, window.size=5)
# tds_5_star 
# 
# tds_5$domain
# tds_5_star$domain
# 
# 



