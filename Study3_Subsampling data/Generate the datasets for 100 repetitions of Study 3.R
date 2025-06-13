
# Repeat the analysis 100 times for Study 3: Subsampling from a bulk Hi-C dataset 

load("C:/Users/yqliu/OneDrive - The Ohio State University/2023Spring/RWR Manuscript/2024Autumn Version/Code upload to GitHub/Datasets/Study3_Subsampling data/Bulk_and_subsampled_data.RData")

source("KR Normalization2.R")

for(iter in 1:100){
  i = 0
  while(exists("subsample_chr22_KR") == FALSE){
    tryCatch({
      subsample_chr22 <- Subsample(bulk_chr22, alpha=0.1, seed=2021+iter*4+i)
      subsample_chr22_KR <- KRnorm2(subsample_chr22)},
      error = function(c) "error"
    )
    i = i+1
  }
  
  save(subsample_chr22, subsample_chr22_KR,
       file=paste0("~/Study3_Subsampling data/Subsampled", iter, " data.RData"))

  rm(subsample_chr22, subsample_chr22_KR)
}
