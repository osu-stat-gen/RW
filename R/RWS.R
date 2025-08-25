
#' Random Walk with limited number of steps (RWS) function
#'
#' Computes the \eqn{k}-step transition matrix \eqn{Q = P^k} on a symmetric
#' doubly-stochastic transition matrix \eqn{P} (here `normmatrix`) by repeated
#' multiplication. This corresponds to a random walk of `step`
#' transitions on the underlying graph/Markov chain.
#'
#' @param normmatrix Numeric symmetric square matrix \eqn{P} with each row/column summing to 1
#'   (doubly-stochastic). The function does not renormalize; ensure the matrix
#'   is normalization beforehand (e.g., by `KRnorm2()` function).
#' @param step Positive integer number of steps \eqn{k} (default: `3`).
#'
#'
#' @return A numeric matrix `Q` of the same dimension as `normmatrix`, equal to \eqn{P^k}.
#'
#' @details
#' This is a straightforward dense computation of \eqn{P^k} via repeated
#' `%*%`. For large matrices and/or large `step`, consider sparse matrices
#' (e.g., via the **Matrix** package) or fast exponentiation.
#'
#' @examples
#' data(sim1.P.KR)
#' # Random walk at 3 steps
#' Q <- RWS(sim1.P.KR, step=3)
#' Q
#'
#' @export
RWS <- function(normmatrix, step = 3) {
  # ---- basic checks ----
  if (!is.matrix(normmatrix)) normmatrix <- as.matrix(normmatrix)
  if (!is.numeric(normmatrix)) stop("'normmatrix' must be numeric.", call. = FALSE)
  if (nrow(normmatrix) != ncol(normmatrix)) stop("'normmatrix' must be square.", call. = FALSE)
  if (length(step) != 1L || !is.finite(step) || step <= 0 || step != as.integer(step)) {
    stop("'step' must be a positive integer.", call. = FALSE)
  }
  if (any(!is.finite(normmatrix))) stop("'normmatrix' contains NA/Inf.", call. = FALSE)

  step <- as.integer(step)

  # check row-stochasticity (warning only)
  rs <- rowSums(normmatrix)
  if (any(abs(rs - 1) > 1e-6)) {
    warning("Some rows of 'normmatrix' do not sum to 1 (row-stochasticity check failed). ",
            "Ensure the matrix is normalized before calling RWS().")
  }

  # Preserve names
  rn <- rownames(normmatrix); cn <- colnames(normmatrix)

  # ---- compute P^step ----
  if (step == 1L) {
    Q <- normmatrix
  } else {
    Q <- normmatrix
    for (k in 2:step) {
      Q <- Q %*% normmatrix
    }
  }

  dimnames(Q) <- list(rn, cn)

  return(Q)
}

