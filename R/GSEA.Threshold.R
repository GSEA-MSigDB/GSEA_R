#' Set thresholds for expression matrix
#'
#' `GSEA.Threshold` sets threshold and ceiling values to pre-process gene expression matrix
#'
#' Internal `GSEA` function.
#'
#' @keywords internal
#'

GSEA.Threshold <-
function(V, thres, ceil) {

 V[V < thres] <- thres
 V[V > ceil] <- ceil
 return(V)
}
