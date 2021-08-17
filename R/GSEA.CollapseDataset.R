#' Maps user supplied identifiers to Gene Symbols
#'
#' `GSEA.CollapseDataset` applies a .CHIP file to the dataset to gene symbols using a specified collapse function
#'
#' Internal `GSEA` function invoked if collapse.dataset == TRUE.
#'
#' @param dataplatform Dataframe returned by `GSEA.ReadCHIPFile`
#' @param gct Dataframe returned by `GSEA.Gct2Frame` after processing to remove unnecessary GCT header fields
#' @param collapse.mode Method for collapsing the dataset, accepts "max", "median", "mean", "sum", (default: NOCOLLAPSE)
#'
#' @return Gene by sample matrix converted into the namespace specified by the applied CHIP file
#'
#' @keywords internal
#'
#' @import dplyr
#' @importFrom rlang .data
#'

GSEA.CollapseDataset <-
function(dataplatform, gct, collapse.mode) {
 probemap <- unique(dataplatform[, c("Probe.Set.ID", "Gene.Symbol")])
 annotate <- unique(dataplatform[, c("Gene.Symbol", "Gene.Title")])
 mappedgct <- merge(x = probemap, y = gct, by.x = 1, by.y = 1, all = FALSE, no.dups = FALSE)
 mappedgct = unique(subset(mappedgct, select = -c(1)))
 mappedgct <- unique(mappedgct[, -c(2)])
 mappedgct[, c(2:ncol(mappedgct))] <- sapply(mappedgct[, c(2:ncol(mappedgct))],
  as.numeric)

 if (collapse.mode == "max")   #MAX
  {
   mappedexp <- mappedgct %>% group_by(.data$Gene.Symbol) %>% summarise_all(max, na.rm = TRUE) %>%
    data.frame()
  }
 if (collapse.mode == "median")  #Median
  {
   mappedexp <- mappedgct %>% group_by(.data$Gene.Symbol) %>% summarise_all(median, na.rm = TRUE) %>%
    data.frame()
  } 
if (collapse.mode == "mean") #Mean
  {
   mappedexp <- mappedgct %>% group_by(.data$Gene.Symbol) %>% summarise_all(mean, na.rm = TRUE) %>%
    data.frame()
  } 
 if (collapse.mode == "sum") #SUM
  {
   mappedexp <- mappedgct %>% group_by(.data$Gene.Symbol) %>% summarise_all(sum, na.rm = TRUE) %>%
    data.frame()
  } 
 mappedexp_2 <- unique(merge(x = annotate, y = mappedexp, by.x = "Gene.Symbol",
  by.y = "Gene.Symbol"))
 colnames(mappedexp_2)[2] <- "Description"
 return(mappedexp_2)
}
