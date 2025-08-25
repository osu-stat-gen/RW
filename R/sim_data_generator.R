
#' Simulate Hi-C count matrices with domain structure and zero inflation
#'
#' Generates a symmetric Hi-C contact matrix with:
#' (1) a distance-decay background following \eqn{\log(\text{count}) \propto \beta_{1} \log(\text{genomic distance})};
#' (2) domain (TAD-like) blocks with extra signal;
#' (3) zero inflation inside and outside domains; and
#' (4) negative-binomial overdispersion.
#'
#' @importFrom stats rnbinom rpois runif median
#'
#' @param n Integer number of bins (matrix dimension). Default `1000`.
#' @param const Positive scaling constant representing diagonal and first
#'   off-diagonal magnitude. Default `20`.
#' @param beta1 Decay exponent in the biophysical law (should be negative),
#'   default `-0.4`.
#' @param disp Positive over-dispersion factor for the Negative Binomial distribution:
#'   we use `size = disp * mu`, so `Var = (1 + 1/disp) * mu`. Default `5`.
#' @param k Integer number of domains (at least 2 and at most `n`). Default `20`.
#' @param a,b Range for uniform domain effect scalars \eqn{\lambda_{b} \sim Unif(a,b)}.
#'   Defaults `a = 0.5`, `b = 1`.
#' @param pi1 Within-domain zero proportion in \[0,1\]. We zero out approximately
#'   `pi1 * 100%` of upper-triangular entries below the domain median, which
#'   yields about `pi1/2` zeros per domain after symmetrization. Default `0.2`.
#' @param pi2 Outside-domain zero proportion in \[0,1\]. We zero out approximately
#'   `pi2 * 100%` of the upper-triangular non-domain entries below their median,
#'   yielding about `pi2/2` zeros outside domains after symmetrization. Default `0.4`.
#'
#' @return A list with:
#' \describe{
#'   \item{simulated.matrix}{Symmetric `n` by `n` integer matrix of simulated counts.}
#'   \item{true.bound}{Integer vector of length `k` giving cumulative domain endpoints
#'         (last element is `n`). The domain intervals are `(0, true.bound[1]]`,
#'         `(true.bound[1], true.bound[2]]`, ...}
#' }
#'
#' @details
#' **Background (distance-decay):** We build a mean matrix `D` with diagonal
#' `const` and the `i`-th off-diagonal proportional to `const * i^{beta1}`.
#' Counts are sampled as Negative Binomial with `mu = D`, `size = disp * D`.
#'
#' **Domain effects:** We draw `k-1` random boundaries and independent
#' `lambda_b ~ U(a,b)` domain scalars, then add Poisson noise with
#' `lambda = lambda_b * domain_block` inside each domain.
#'
#' **Zero inflation:** We zero a fraction of small (below-median) entries in
#' the **upper triangle** inside each domain (rate `pi1`) and in the outside
#' non-domain region (rate `pi2`), then mirror to enforce symmetry.
#'
#' #' @examples
#' set.seed(123)
#' sim1 <- sim_data_generator(n = 1000, const = 20, beta1 = -0.4, disp = 5, k = 20,
#'                            a = 0.5, b = 1, pi1 = 0.2, pi2 = 0.4)
#' dim(sim1$simulated.matrix)
#'
#' @export
sim_data_generator <- function(n = 1000, const = 20, beta1 = -0.4, disp = 5, k = 20,
                               a = 0.5, b = 1, pi1 = 0.2, pi2 = 0.4){

  # Simulate a matrix following the biophysical law:
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
