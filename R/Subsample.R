
#' Subsample a contact matrix by multinomial thinning
#'
#' Randomly subsample approximately \eqn{\alpha \times 100\%} of the total
#' contacts from a symmetric Hi-C contact matrix. The subsample is drawn
#' proportionally to the original counts over the upper triangle (including
#' the main diagonal), then mapped back to a symmetric matrix.
#'
#' @importFrom stats rmultinom
#'
#' @param A Numeric, square, symmetric Hi-C contact matrix (non-negative counts).
#' @param alpha Fraction in \eqn{(0, 1]} giving the target proportion of total
#'   contacts to retain (default `0.1`).
#' @param seed Optional integer seed for reproducibility (default `2022`).
#'
#' @return A symmetric matrix of the same dimension as `A`, containing the
#'   subsampled counts. In expectation, each cell is scaled by `alpha`.
#'
#' @details
#' Let \eqn{N=\sum_{i \le j} A_{ij}} be the total contacts in the upper triangle
#' (including diagonal). The procedure draws a single multinomial sample of size
#' \eqn{\lceil \alpha N \rceil} with cell probabilities
#' \eqn{p_{ij} = A_{ij}/N} over the upper triangle, then mirrors to enforce
#' symmetry. This is equivalent to *multinomial thinning*, so
#' \eqn{\mathbb{E}[\tilde A_{ij}] \approx \alpha A_{ij}}.
#'
#' @examples
#' data(sim3.bulk)
#' data(sim3.sub)
#' sub <- Subsample(sim3.bulk, alpha=0.1, seed=2022)
#' # Check if the subsampled data matches with the one produced using the same seed
#' sum(abs(sub - sim3.sub))    # should = 0
#'
#' @export
Subsample <- function(A, alpha = 0.1, seed = 2022) {
  ## ---- checks -------------------------------------------------------------
  if (!is.matrix(A)) A <- as.matrix(A)
  if (!is.numeric(A)) stop("'A' must be numeric.", call. = FALSE)
  if (nrow(A) != ncol(A)) stop("'A' must be square.", call. = FALSE)
  if (any(A < 0, na.rm = TRUE)) stop("'A' must be nonnegative.", call. = FALSE)
  if (!is.finite(alpha) || alpha <= 0 || alpha > 1)
    stop("'alpha' must be in (0, 1].", call. = FALSE)

  n <- nrow(A)
  rn <- rownames(A); cn <- colnames(A)

  ## Total sequencing depth from upper triangle (diagonal included)
  totaldepth <- sum(A[upper.tri(A, diag = TRUE)], na.rm = TRUE)

  # Trivial cases
  if (totaldepth == 0 || alpha == 0) {
    out <- matrix(0, n, n)
    dimnames(out) <- list(rn, cn)
    return(out)
  }

  ## ---- Multinomial thinning ----------------------------------------------
  probs <- matrixtovec(A) / totaldepth
  # guard against numerical drift (e.g., all zeros or rounding)
  if (!all(is.finite(probs)) || sum(probs) <= 0)
    stop("Invalid probability vector derived from 'A'. Check that 'A' has positive mass.", call. = FALSE)

  if (!is.null(seed)) {
    set.seed(as.integer(seed))
  } else { set.seed(2022) }

  size <- as.integer(ceiling(alpha * totaldepth))
  sam <- as.vector(stats::rmultinom(n = 1, size = size, prob = probs))

  ## ---- Reconstruct symmetric matrix --------------------------------------
  sub <- vectomatrix(sam)

  # Keep dimension names
  dimnames(sub) <- list(rn, cn)

  return(sub)
}

