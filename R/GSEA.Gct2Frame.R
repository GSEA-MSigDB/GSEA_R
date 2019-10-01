#' Reads .GCT file into a data frame for processing
#'
#' `GSEA.Gct2Frame` is a wrapper for `read.table` with specified parameters for reading in a .GCT file
#'
#' Internal `GSEA` function.
#'
#' @keywords internal
#'

GSEA.Gct2Frame <- function(filename = "NULL") {
 ds <- read.table(filename, sep = "\t", comment.char = "", quote = "", stringsAsFactors = FALSE, 
  fill = TRUE, header = F)
 ds <- ds[-c(1), ]
 ds <- ds[-c(1), ]
 colnames(ds) <- ds[c(1), ]
 ds <- ds[-c(1), ]
 return(ds)
}
