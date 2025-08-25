
#' Random Walk with Restart (RWR) function
#'
#' Computes the steady-state random-walk-with-restart matrix \eqn{Q} on a symmetric
#' doubly-stochastic transition matrix \eqn{P} (here `normmatrix`), by fixed-point
#' iteration:
#' \deqn{Q_{t+1} = (1-\alpha)\, Q_t P + \alpha I, \qquad Q_0 = I.}
#' At convergence, \eqn{Q} solves \eqn{Q = \alpha(I - (1-\alpha)P)^{-1}} (when
#' the inverse exists). Each column of \eqn{Q} corresponds to the diffusion
#' starting at one node with restart probability \eqn{\alpha}.
#'
#' @param normmatrix Numeric symmetric square matrix \eqn{P} with each row/column summing to 1
#'   (doubly-stochastic). This function does not renormalize; ensure the matrix
#'   is normalized beforehand (e.g., by `KRnorm2()` function).
#' @param alpha Restart probability in \eqn{(0,1)}. Default: `0.1`.
#' @param max_iter Maximum number of iterations. Default: `300`.
#' @param tol Convergence tolerance on the Frobenius norm of the iterate
#'   difference. Default: `1e-6`.
#'
#' @return A numeric matrix `Q` of the same dimension as `normmatrix`.
#'
#' @details
#' For better numerical behavior, ensure `normmatrix` is nonnegative, symmetric, and
#' stochastic (`rowSums(normmatrix) == 1` within small tolerance). The
#' iteration stops when `||Q_{t+1} - Q_t||_F < tol` or when `max_iter` is reached.
#'
#' @examples
#' data(sim1.P.KR)
#' # Random walk with restart where restart probability is 0.1
#' Q <- RWR(sim1.P.KR, alpha=0.1, max_iter=300, tol=1e-6)
#'
#' @export
RWR <- function(normmatrix, alpha = 0.1, max_iter = 300, tol = 1e-6) {
  # ---- basic checks ----
  if (!is.matrix(normmatrix)) normmatrix <- as.matrix(normmatrix)
  if (!is.numeric(normmatrix)) stop("'normmatrix' must be numeric.", call. = FALSE)
  if (nrow(normmatrix) != ncol(normmatrix)) stop("'normmatrix' must be square.", call. = FALSE)
  if (!is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("Restart probability 'alpha' must be strictly between 0 and 1.", call. = FALSE)
  }
  if (max_iter < 1L) stop("'max_iter' must be >= 1.", call. = FALSE)
  if (tol <= 0) stop("'tol' must be > 0.", call. = FALSE)
  if (any(!is.finite(normmatrix))) stop("'normmatrix' contains NA/Inf/NaN.", call. = FALSE)

  # check row-stochasticity (warning only)
  n <- nrow(normmatrix)
  rs <- rowSums(normmatrix)
  if (any(abs(rs - 1) > 1e-6)) {
    warning("Some rows of 'normmatrix' do not sum to 1 (row-stochasticity check failed). ",
            "Ensure you've normalized the matrix before calling RWR().")
  }

  # Preserve names
  rn <- rownames(normmatrix); cn <- colnames(normmatrix)

  # ---- fixed-point iteration ----
  Q <- diag(n)
  delta <- Inf
  for (i in seq_len(max_iter)) {
    Q_new <- (1 - alpha) * (Q %*% normmatrix) + alpha * diag(n)
    delta <- base::norm(Q - Q_new, type = "F")
    Q <- Q_new
    if (delta < tol) break
  }

  message(sprintf("Iterations: %d; last change (Frobenius norm): %.3e", i, delta))
  dimnames(Q) <- list(rn, cn)

  return(Q)
}

