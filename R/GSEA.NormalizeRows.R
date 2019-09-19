#' Standardize rows of a gene expression matrix
#'
#' `GSEA.NormalizeRows` performs row-wise adjustment of a gene expression matrix using the col.mean and standard deviation
#'
#' Internal `GSEA` function.
#'
#' @keywords internal
#'

GSEA.NormalizeRows <-
function(V) {

 row.mean <- apply(V, MARGIN = 1, FUN = mean)
 row.sd <- apply(V, MARGIN = 1, FUN = sd)
 row.n <- length(V[, 1])
 for (i in 1:row.n) {
  if (row.sd[i] == 0) {
   V[i, ] <- 0
  } else {
   V[i, ] <- (V[i, ] - row.mean[i])/row.sd[i]
  }
 }
 return(V)
}
