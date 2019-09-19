#' Read phenotype assignments from .CLS file
#'
#' `GSEA.ReadClsFile` reads a class vector CLS file to define dataset relationships 
#'
#' Read in a class vector CLS file to defines phenotype and class label
#' vectors for the samples in a gene expression file (GCT format).
#'
#' @keywords internal
#'

GSEA.ReadClsFile <-
function(file = "NULL") {

 cls.cont <- readLines(file)
 num.lines <- length(cls.cont)
 cls.cont[[3]] <- gsub("\\t", " ", cls.cont[[3]])  #Converts any tabs to spaces
 class.list <- unlist(strsplit(cls.cont[[3]], " "))  #Splits CLS on spaces
 s <- length(class.list)
 t <- table(class.list)
 l <- length(t)
 phen <- vector(length = l, mode = "character")
 phen.label <- vector(length = l, mode = "numeric")
 class.v <- vector(length = s, mode = "numeric")
 for (i in 1:l) {
  phen[i] <- noquote(names(t)[i])
  phen.label[i] <- i - 1
 }
 for (i in 1:s) {
  for (j in 1:l) {
   if (class.list[i] == phen[j]) {
    class.v[i] <- phen.label[j]
   }
  }
 }
 return(list(phen = phen, class.v = class.v))
}
