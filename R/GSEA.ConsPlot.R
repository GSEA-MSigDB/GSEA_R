#' Plots a heatmap of a consensus matrix
#'
#' `GSEA.ConsPlot` plots heatmaps of the consensus matrix from leading edge analysis
#'
#' Invoked by `GSEA.Analyze.Sets` to plot the consense matrix of the leading edge analysis.
#'
#' @keywords internal
#'

GSEA.ConsPlot <-
function(V, col.names, main = " ", sub = " ", xlab = " ", ylab = " ") {

 cols <- length(V[1, ])
 B <- matrix(0, nrow = cols, ncol = cols)
 max.val <- max(V)
 min.val <- min(V)
 for (i in 1:cols) {
  for (j in 1:cols) {
   k <- cols - i + 1
   B[k, j] <- max.val - V[i, j] + min.val
  }
 }



 # col.map <- c(rainbow(100, s = 1.0, v = 0.75, start = 0.0, end = 0.75, gamma =
 # 1.5), '#BBBBBB', '#333333', '#FFFFFF')
 col.map <- rev(c("#0000FF", "#4040FF", "#7070FF", "#8888FF", "#A9A9FF", "#D5D5FF",
  "#EEE5EE", "#FFAADA", "#FF9DB0", "#FF7080", "#FF5A5A", "#FF4040", "#FF0D1D"))

 # max.size <- max(nchar(col.names))
 par(mar = c(5, 15, 15, 5))
 image(1:cols, 1:cols, t(B), col = col.map, axes = FALSE, main = main, sub = sub,
  xlab = xlab, ylab = ylab)

 for (i in 1:cols) {
  col.names[i] <- substr(col.names[i], 1, 25)
 }
 col.names2 <- rev(col.names)

 size.col.char <- ifelse(cols < 15, 1, sqrt(15/cols))

 axis(2, at = 1:cols, labels = col.names2, adj = 0.5, tick = FALSE, las = 1, cex.axis = size.col.char,
  font.axis = 1, line = -1)
 axis(3, at = 1:cols, labels = col.names, adj = 1, tick = FALSE, las = 3, cex.axis = size.col.char,
  font.axis = 1, line = -1)

 return()
}
