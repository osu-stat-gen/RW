
#' Heatmap visualization of Hi-C data with optional domain boundaries overlay
#'
#' Plot a symmetric Hi-C contact matrix as a heatmap (upper-left origin),
#' with optional truncation of large values and optional overlay of
#' domain boundaries.
#'
#' @import ggplot2
#' @importFrom grid unit
#' @importFrom rlang .data
#'
#' @param hic Numeric, square matrix of Hi-C contact counts (non-negative).
#' @param trunc Logical; if `TRUE`, values larger than `thres` are clipped
#'   before plotting. Default `FALSE`.
#' @param thres Positive number used when `trunc = TRUE` to cap large values.
#'   Default `20`.
#' @param title Character plot title. Default empty string.
#' @param add.bound Logical; if `TRUE`, draw boundary line segments that mark
#'   blocks defined by `bound`. Default `FALSE`.
#' @param bound Optional increasing integer vector of domain boundary positions
#'   (end indices of consecutive blocks). For example, if `nrow(hic) = 100`
#'   and `bound = c(20, 55, 100)`, the diagonal blocks are `(1,...,20)`, `(21,...,55)`,
#'   `(56,...,100)`. If `add.bound = TRUE`, `bound` must be supplied.
#' @param col Color for boundary lines. Default `"blue"`.
#'
#' @return A `ggplot` object (heatmap), optionally with boundaries overlaid.
#'
#' @details
#' For display, the matrix is transposed and the y-axis is flipped so that
#' bin `(1,1)` appears in the top-left (conventional Hi-C view).
#'
#' If `bound` vector does not include the last bin (`nrow(hic)`), it will
#' be appended internally for drawing the complete outline.
#'
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   data(sim2.E)
#'   data(sim2.true.bound.E)
#'
#'   # Heatmap without TAD boundaries
#'   p1 <- hic_heatmap(sim2.E, trunc=TRUE, thres=30, add.bound=FALSE, title="Sim2 heatmap")
#'
#'   # Heatmap with true TAD boundaries
#'   p2 <- hic_heatmap(sim2.E, trunc=TRUE, thres=30, title="Sim2 heatmap with TADs",
#'                     add.bound=TRUE, bound = sim2.true.bound.E, col="blue")
#' }
#'
#' @export
hic_heatmap <- function(hic, trunc = FALSE, thres = 20, title = "",
                        add.bound = FALSE, bound = NULL, col = "blue") {

  ## ---- basic checks -------------------------------------------------------
  if (!is.matrix(hic)) hic <- as.matrix(hic)
  if (!is.numeric(hic)) stop("'hic' must be numeric.", call. = FALSE)
  if (nrow(hic) != ncol(hic)) stop("'hic' must be square.", call. = FALSE)
  if (!is.logical(trunc) || length(trunc) != 1L) stop("'trunc' must be TRUE/FALSE.", call. = FALSE)
  if (!is.finite(thres) || thres <= 0) stop("'thres' must be a positive number.", call. = FALSE)
  if (!is.logical(add.bound) || length(add.bound) != 1L) stop("'add.bound' must be TRUE/FALSE.", call. = FALSE)
  n <- nrow(hic)

  ## ---- data prep (flip Y so origin is top-left) ---------------------------
  x <- seq_len(n); y <- seq_len(n)
  df <- expand.grid(X = x, Y = y)
  vals <- as.vector(t(hic)[, n:1, drop = FALSE])
  if (trunc) vals <- pmin(vals, thres)
  df$bulk <- vals

  ## ---- base heatmap -------------------------------------------------------
  bulkplot <- ggplot(df, aes(.data$X, .data$Y, fill = .data$bulk)) +
    geom_tile() +
    theme_bw() +
    scale_fill_gradient(low = "white", high = "red3",
                        limits = c(0, max(df$bulk, na.rm = TRUE)), name = "") +
    theme(panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank()) +
    labs(x = "", y = "", title = title) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.key.size = unit(1.5, "cm"),
          legend.key.height = unit(1.5, "cm"),
          legend.key.width = unit(0.8, "cm"),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 14),
          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) +
    coord_equal()

  ## ---- add boundaries (optional) -----------------------------------------
  if (!add.bound) {
    return(bulkplot)
  } else {
    if (is.null(bound) || length(bound) < 1L)
      stop("'bound' must be provided when add.bound = TRUE.", call. = FALSE)
    if (any(!is.finite(bound)) || any(bound < 1) || any(bound > n))
      stop("'bound' must contain valid indices in 1..nrow(hic).", call. = FALSE)
    bound <- sort(unique(as.integer(bound)))
    # ensure final boundary at n
    if (max(bound) < n) bound <- c(bound, n)

    return(add_bound(bulkplot, bound = bound, col = col))
  }
}


