
# Repeat the analysis 100 times for Study 2: Biophysical-law-based simulated Hi-C data 

source("KR Normalization2.R")

for(iter in 1:100){
  set.seed(2022416+iter)
  n = 1000
  sim1 <- sim_data_generator(n = 1000, pi1 = 0.1, pi2 = 0.4)

  # The simulated matrix
  E <- sim1$simulated.matrix
  # The true domain boundary endpoints
  true.bound.E <- sim1$true.bound

  E_KR <- KRnorm2(E)

  save(E, E_KR, true.bound.E,
       file=paste0("~/Study2_Biophysical-law data/Iter", iter, " simulated data.RData"))
}
