#' Produces overlap heatmap for leading edge analysis
#'
#' `GSEA.HeatMapPlot2` plots a heatmap of a gene set leading edge overlap matrix
#'
#' Internal function invoked by `GSEA.Analyze.Sets` to plot heatmaps of leading edge overlaps.
#'
#' @keywords internal
#'

GSEA.HeatMapPlot2 <-
function(V, row.names = "NA", col.names = "NA", main = " ",
 sub = " ", xlab = " ", ylab = " ", color.map = "default") {

 n.rows <- length(V[, 1])
 n.cols <- length(V[1, ])

 if (color.map == "default") {
  color.map <- rev(rainbow(100, s = 1, v = 0.75, start = 0, end = 0.75))
 }

 heatm <- matrix(0, nrow = n.rows, ncol = n.cols)
 heatm[1:n.rows, ] <- V[seq(n.rows, 1, -1), ]

 par(mar = c(7, 15, 5, 5))
 image(1:n.cols, 1:n.rows, t(heatm), col = color.map, axes = FALSE, main = main,
  sub = sub, xlab = xlab, ylab = ylab)

 if (length(row.names) > 1) {
  size.row.char <- ifelse(n.rows < 15, 1, sqrt(15/n.rows))
  size.col.char <- ifelse(n.cols < 15, 1, sqrt(25/n.cols))
  # size.col.char <- ifelse(n.cols < 2.5, 1, sqrt(2.5/n.cols))
  for (i in 1:n.rows) {
   row.names[i] <- substr(row.names[i], 1, 40)
  }
  row.names <- row.names[seq(n.rows, 1, -1)]
  axis(2, at = 1:n.rows, labels = row.names, adj = 0.5, tick = FALSE, las = 1,
   cex.axis = size.row.char, font.axis = 1, line = -1)
 }

 if (length(col.names) > 1) {
  axis(1, at = 1:n.cols, labels = col.names, tick = FALSE, las = 3, cex.axis = size.col.char,
   font.axis = 2, line = -1)
 }
 return()
}
