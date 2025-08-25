
#' Knight–Ruiz (KR) normalization with zero-row handling
#'
#' Performs matrix balancing via the Knight–Ruiz algorithm (KR normalization),
#' adapted for contact matrices that may contain all-zero rows/columns and
#' missing values. Zero rows/columns are temporarily removed for the solve and
#' inserted back afterward with their diagonal entries set to 1 (others 0).
#'
#' @param B A numeric, square (n by n) matrix (e.g., a Hi-C contact map). May
#'   contain \code{NA}s and rows/columns that sum to zero (after treating
#'   \code{NA} as 0 for the purpose of zero-detection).
#'
#' @return A numeric matrix of the same dimension as \code{B}, KR-normalized.
#'   Any zero rows/columns are restored with identity on their diagonal.
#'   \code{NA} entries in the nonzero block are reintroduced in the output.
#'
#' @details
#' The core iterations follow Knight & Ruiz (2012). Internally, the method
#' solves for a positive scaling vector \eqn{x} such that
#' \eqn{\mathrm{diag}(x) \, B \, \mathrm{diag}(x)} is doubly
#' stochastic (on the submatrix excluding zero rows/columns). We treat
#' \code{NA}s as zeros during the solve but place them back into the result
#' after scaling.
#'
#' @references
#' Knight, P. A., & Ruiz, D. (2012). A fast algorithm for matrix balancing.
#' \emph{IMA Journal of Numerical Analysis}, 32(3), 906–929.
#' \doi{10.1093/imanum/drs019}
#'
#' @examples
#' set.seed(123)
#' M <- matrix(rpois(25, 5), 5, 5)
#' M[lower.tri(M)] <- t(M)[lower.tri(M)]
#' diag(M) <- 0
#' M[4, ] <- 0; M[, 4] <- 0              # a zero row/column
#' M.KR <- KRnorm2(M)
#' M.KR
#'
#' @export
KRnorm2 <- function(B) {
  # create a placeholder matrix
  C = B

  # remove any cols/rows of 0s
  zeros = unique(which(colSums(B) == 0), which(rowSums(B) == 0))

  # save col/row names
  cnames <- colnames(B)

  if (length(zeros) > 0) {
    A = B[-zeros, -zeros]
    message(paste0('Cols/Rows being zeros: ', length(zeros), ' of them'))
    message(paste(" ", zeros, sep = " "))
  } else {
    A = B
  }

  # initialize
  tol = 1e-6; delta = 0.1; Delta = 3; fl = 0;
  # change NAs in matrix to 0's
  NAlist = which(is.na(A), arr.ind = TRUE)
  A[is.na(A)] = 0
  n = nrow(A)
  e = matrix(1, nrow=n, ncol = 1)
  x0 = e
  res = matrix(nrow = n, ncol=1)
  # inner stopping criterior
  g=0.9; etamax = 0.1;
  eta = etamax; stop_tol = tol*.5;
  x = x0; rt = tol^2; v = x*(A %*% x); rk = 1-v;
  rho_km1 = t(rk) %*% rk; rout = rho_km1; rold = rout;

  MVP = 0; # We'll count matrix vector products.
  i = 0; # Outer iteration count.

  while(rout > rt) { # Outer iteration
    i = i + 1; k = 0; y = e;
    innertol = max(c(eta^2 %*% rout, rt));

    while( rho_km1 > innertol ) { #Inner iteration by CG
      k = k+1
      if(k==1) {
        z = rk/v; p = z; rho_km1 = t(rk)%*%z;
      }else {
        beta = rho_km1 %*% solve(rho_km2)
        p = z + beta*p
      }

      # update search direction efficiently
      w = x * (A%*%(x*p)) + v*p
      alpha = rho_km1 %*% solve(t(p) %*% w)
      ap = c(alpha) * p

      # test distance to boundary of cone
      ynew = y + ap;
      if(min(ynew) <= delta) {
        if(delta == 0) break()
        ind = which(ap < 0);
        gamma = min((delta - y[ind])/ap[ind]);
        y = y + gamma * ap;
        break()
      }
      if(max(ynew) >= Delta) {
        ind = which(ynew > Delta);
        gamma = min((Delta-y[ind])/ap[ind]);
        y = y + gamma * ap;
        break()
      }
      y = ynew;
      rk = rk - c(alpha) * w; rho_km2 = rho_km1;
      Z = rk/v; rho_km1 = t(rk) %*% z;
    }

    x = x*y; v = x*(A %*% x);
    rk = 1 - v; rho_km1 = t(rk) %*% rk; rout = rho_km1;
    MVP = MVP + k + 1;

    # Update inner iteration stopping criterion.
    rat = rout %*% solve(rold); rold = rout; res_norm = sqrt(rout);
    eta_o = eta; eta = g %*% rat;
    if(g %*% eta_o^2 > 0.1) {
      eta = max(c(eta, g %*% eta_o^2));
    }
    eta = max(c(min(c(eta, etamax)), stop_tol/res_norm));
  }

  result = t(t(x[,1]*A)*x[,1])

  # reintroduce NAs in final matrix
  if(nrow(NAlist) > 0) {
    idx <- as.matrix(NAlist[, 1:2])
    result[idx] <- NA
  }

  # add zero rows/columns back
  if (length(zeros) > 0) {
    C[-zeros, -zeros] = result
    C[zeros, zeros] = diag(1, nrow=length(zeros), ncol=length(zeros))
  } else {
    C = result
  }

  return(C)
}

