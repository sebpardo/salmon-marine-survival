
# Non-centered parameterization comparison
compare_cp_ncp <- function(cp_plot, ncp_plot, ncol = 2, ...) {
  bayesplot_grid(
    cp_plot, ncp_plot,
    grid_args = list(ncol = ncol),
    subtitles = c("Centered parameterization",
                  "Non-centered parameterization"),
    ...
  )
}

# function from http://www.sthda.com/english/wiki/visualize-correlation-matrix-using-correlogram
# to estimate p-values for all correlations in matrix
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}