library(matrixcalc)

# A function that generates a simulated Hi-C data with 
# 1. given number of bins (n), 
# 2. a scaling constant (const) which decides the maximum count in the diagonal and first off-diagonal elements in the Hi-C matrix, 
# 3. the coefficient in the biophysical law (beta1, which must be negative), 
# 4. level of over-dispersion in the negative binomial distribution (disp, theta = disp*mu, so var = (1+1/disp)*mu), 
# 5. number of domains (k, at least 2 and at most n), 
# 6. the range of lambda (min = a, max = b), where lambda is a vector of length k (sampled from Uniform(a,b)) which adds the domain effect to the background matrix, 
# 7. Within-domain zero proportion (pi1): pi1*100% of the counts that are inside domain b and less than the median count of domain b are set to be zeros (so each domain has pi1/2 of zeros), b = 1,...,k,  
# 8. Outside-domain zero proportion (pi2): pi2*100% of the counts that are less than the median count for all the regions that are not inside any domain are set to be zeros (so the outside domain region has pi2/2 of zeros). 


sim_data_generator <- function(n = 1000, const = 20, beta1 = -0.4, disp = 5, k = 20,
                               a = 0.2, b = 0.5, pi1 = 0.1, pi2 = 0.4){
  
  # Simulate a matrix with the biophysical law:
  # log(count) \propto (beta1) log(genomic distance)
  # This will be a mean matrix used to generate the counts 
  D <- matrix(0, nrow=n, ncol=n)
  diag(D) <- rep(const/2, n)
  
  for(i in 1:(n-1)){
    D[seq(1+n*i, n*n, 1+n)] <- const*i^(beta1)
  }
  D <- (D+t(D))
  
  
  # Generate negative binomial distributed random counts with means as elements in D
  # for example, when theta = 5*mu, var = 1.2*mean 
  # This is a background matrix -- the domain effect is added later
  E <- matrix(rnbinom(n=n^2, mu=D, size=disp*D), nrow=n, ncol=n)
  E[lower.tri(E, diag=FALSE)] <- 0
  diag(E) <- diag(E)/2
  E <- (E+t(E))
  
  # Generate k-1 random domain boundary points, the last domain ends at n   
  true.bound.E <- c(sort(sample(1:(n-1), size=(k-1), replace=FALSE)), n)

  # Generate k random scaled domain effect sizes  
  lambda.b.E <- runif(k, min=a, max=b)
  
  # Add Poisson random counts to each of the domains and introduce zero counts (with proportion pi1)
  bound <- c(0, true.bound.E)
  for(i in 1:k){
    # find the i^th domain and its size 
    domain.i <- E[(bound[i]+1):bound[i+1], (bound[i]+1):bound[i+1]]
    domain.i.size <- (bound[i+1]-bound[i])
    
    # Enhance the i^th domain 
    noise.i <- matrix(rpois(n=domain.i.size^2, lambda=lambda.b.E[i]*domain.i), nrow=domain.i.size)
    domain.i <- domain.i + noise.i
    
    # sample pi1*100% of the counts that are less than the median to be zeros
    candidate.index <- which(domain.i < median(domain.i))
    candidate.index <- candidate.index[which( ((candidate.index+domain.i.size-1)%%domain.i.size)+1 <= (candidate.index+domain.i.size-1)%/%domain.i.size )]  # row index <= column index, meaning in the upper triangle
    zero.index <- sample(candidate.index, size=floor(pi1*0.5*domain.i.size*(domain.i.size+1)), replace=FALSE)
    domain.i[zero.index] <- 0
    
    E[(bound[i]+1):bound[i+1], (bound[i]+1):bound[i+1]] <- domain.i
    
    # inside a domain only the upper triangle is completed for now,
    # later we will match the lower triangle of the whole count matrix to the upper triangle 
    # make the counts inside a domain symmetric
  }
  
  
  # Introduce zero counts to the non-domain regions (with proportion pi2)
  row.index <- c()
  col.index <- c()
  for(i in 2:(length(bound)-1)){
    row.index <- c(row.index, rep(1:bound[i], times=(bound[i+1]-bound[i])))
    col.index <- c(col.index, rep((bound[i]+1):bound[i+1], each=bound[i]))
  }
  candidate.index <- ((row.index-1)*n+col.index)
  candidate.index <- candidate.index[which(E[candidate.index] < median(E[candidate.index]))]
  zero.index <- sample(candidate.index, size=floor(pi2*0.5*length(row.index)), replace=FALSE)
  E[zero.index] <- 0 
  
  # make the whole count matrix symmetric 
  E[lower.tri(E, diag=FALSE)] <- 0
  diag(E) <- diag(E)/2
  E <- (E+t(E))
  
  result <- list(simulated.matrix = E,
                 true.bound = true.bound.E)
  
  return(result)
}










# ##  An Example  ##
# 
# # Set seed to make the result repeatable 
# set.seed(2022417)
# # Use the default parameter settings
# result <- sim_data_generator(n=300)
# # The simulated matrix
# E <- result$simulated.matrix
# # The true domain boundary endpoints
# true.bound.E <- result$true.bound
# 
# 
# # Visualization # 
# 
# # Function that takes the heatmap and boundary points as input, and returns the heatmap with lines marking the domains #
# add_bound <- function(normplot, bound, col="red"){
#   plot1 <- normplot
#   bound <- c(0, bound[bound>0])
#   nb <- length(bound)
#   if(nb<=2){
#     segment = data.frame(
#       x = c(rep(0.5,3), bound[2]+0.5),
#       y = c(bound[2]+0.5, 0.5, rep(bound[2]+0.5, 2)), 
#       xend = c(rep(bound[2]+0.5, 2), 0.5, bound[2]+0.5), 
#       yend = c(bound[2]+0.5, rep(0.5,3))
#     )
#     plot1 <- plot1 + annotate("segment", x = segment$x, xend = segment$xend, y = segment$y, yend = segment$yend, color=col, size=0.1)
#     return(plot1)
#   } else{
#     segment = data.frame(
#       x = c(bound[2:(nb-1)]+0.5, bound[1:(nb-2)]+0.5),
#       y = c(bound[nb]+0.5-bound[1:(nb-2)], bound[nb]+0.5-bound[2:(nb-1)]), 
#       xend = c(bound[2:(nb-1)]+0.5, bound[3:nb]+0.5),
#       yend = c(bound[nb]+0.5-bound[3:nb], bound[nb]+0.5-bound[2:(nb-1)])
#     )
#     plot1 <- plot1 + annotate("segment", x = segment$x, xend = segment$xend, y = segment$y, yend = segment$yend, color=col, size=0.1)
#     return(plot1)
#   }
# }
# 
# 
# # Plot the simulated count matrix #
# library(ggplot2)
# library(cowplot)
# 
# n <- nrow(E)
# x <- 1:n
# y <- 1:n
# data <- expand.grid(X=x, Y=y)
# data$E <- as.vector(t(E)[,n:1])
# 
# Eplot <- ggplot(data, aes(X, Y, fill=E)) + geom_tile() + theme_bw()+
#   scale_fill_gradient(low="white", high="darkblue") +
#   theme(panel.grid=element_blank(), panel.border=element_blank()) +
#   theme(axis.text = element_blank()) + 
#   theme(axis.ticks = element_blank()) +
#   labs(x = "", y = "", title = "Contact count matrix E") +
#   theme(plot.title = element_text(hjust = 0.5),
#         legend.key.size = unit(1, 'cm'), 
#         legend.key.height = unit(0.6, 'cm'), 
#         legend.key.width = unit(0.4, 'cm'), 
#         legend.title = element_text(size=10), 
#         legend.text = element_text(size=8)) +
#   coord_equal() 
# 
# Eplot_true <- Eplot + labs(title = "E + true domains") 
# 
# plot_grid(Eplot, add_bound(Eplot_true, true.bound.E, col="red"), nrow=1)
# 
# 
# 
# 
# # Plot the count vs distance # 
# 
# row.idx <- c()  # index for locus i 
# col.idx <- c()  # index for locus j
# for(i in 1:n){
#   row.idx <- c(row.idx, 1:i)
#   col.idx <- c(col.idx, rep(i, i))
# }
# 
# plot(abs(row.idx - col.idx), E[upper.tri(E, diag=TRUE)], xlab="distance", ylab="count", main="")
# 
# plot(abs(n+1 - data$Y - data$X), data$E, xlab="distance", ylab="count", main="")
