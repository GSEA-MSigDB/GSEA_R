#' Variation filter pre-processing for gene expression matrix
#'
#' `GSEA.VarFilter` apply a variation filter to expression matrix
#'
#' Unused function.
#'
#' @keywords internal
#'

GSEA.VarFilter <-
function(V, fold, delta, gene.names = "NULL") {

 cols <- length(V[1, ])
 rows <- length(V[, 1])
 row.max <- apply(V, MARGIN = 1, FUN = max)
 row.min <- apply(V, MARGIN = 1, FUN = min)
 flag <- array(dim = rows)
 flag <- (row.max/row.min > fold) & (row.max - row.min > delta)
 size <- sum(flag)
 B <- matrix(0, nrow = size, ncol = cols)
 j <- 1
 if (gene.names == "NULL") {
  for (i in 1:rows) {
   if (flag[i]) {
    B[j, ] <- V[i, ]
    j <- j + 1
   }
  }
  return(B)
 } else {
  new.list <- vector(mode = "character", length = size)
  for (i in 1:rows) {
   if (flag[i]) {
    B[j, ] <- V[i, ]
    new.list[j] <- gene.names[i]
    j <- j + 1
   }
  }
  return(list(V = B, new.list = new.list))
 }
}
