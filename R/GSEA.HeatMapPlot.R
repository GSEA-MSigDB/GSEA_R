#' Produce Heat Map for Genes in Dataset
#'
#' `GSEA.HeatMapPlot` plots a heatmap 'pinkogram' of a gene expression matrix including phenotype vector and gene, sample and phenotype labels
#'
#' Internal `GSEA` function invoked if gsea.type="GSEA"
#'
#' @keywords internal
#'

GSEA.HeatMapPlot <-
function(V, row.names = F, col.labels, col.classes, col.names = F,
 main = " ", xlab = " ", ylab = " ") {

 n.rows <- length(V[, 1])
 n.cols <- length(V[1, ])
 row.mean <- apply(V, MARGIN = 1, FUN = mean)
 row.sd <- apply(V, MARGIN = 1, FUN = sd)
 row.n <- length(V[, 1])
 for (i in 1:n.rows) {
  if (row.sd[i] == 0) {
   V[i, ] <- 0
  } else {
   V[i, ] <- (V[i, ] - row.mean[i])/(0.5 * row.sd[i])
  }
  V[i, ] <- ifelse(V[i, ] < -6, -6, V[i, ])
  V[i, ] <- ifelse(V[i, ] > 6, 6, V[i, ])
 }

 mycol <- c("#0000FF", "#0000FF", "#4040FF", "#7070FF", "#8888FF", "#A9A9FF",
  "#D5D5FF", "#EEE5EE", "#FFAADA", "#FF9DB0", "#FF7080", "#FF5A5A", "#FF4040",
  "#FF0D1D", "#FF0000")  # blue-pinkogram colors. The first and last are the colors to indicate the class vector (phenotype). This is the 1998-vintage, pre-gene cluster, original pinkogram color map

 mid.range.V <- mean(range(V)) - 0.1
 heatm <- matrix(0, nrow = n.rows + 1, ncol = n.cols)
 heatm[1:n.rows, ] <- V[seq(n.rows, 1, -1), ]
 heatm[n.rows + 1, ] <- ifelse(col.labels == 0, 7, -7)
 image(1:n.cols, 1:(n.rows + 1), t(heatm), col = mycol, axes = FALSE, main = main,
  xlab = xlab, ylab = ylab)

 if (length(row.names) > 1) {
  numC <- nchar(row.names)
  size.row.char <- 25/(n.rows + 5)
  size.col.char <- 25/(n.cols + 5)
  maxl <- floor(n.rows/1.6)
  for (i in 1:n.rows) {
   row.names[i] <- substr(row.names[i], 1, maxl)
  }
  row.names <- c(row.names[seq(n.rows, 1, -1)], "Class")
  axis(2, at = 1:(n.rows + 1), labels = row.names, adj = 0.5, tick = FALSE,
   las = 1, cex.axis = size.row.char, font.axis = 2, line = -1)
 }

 if (length(col.names) > 1) {
  axis(1, at = 1:n.cols, labels = col.names, tick = FALSE, las = 3, cex.axis = size.col.char,
   font.axis = 2, line = -1)
 }

 C <- split(col.labels, col.labels)
 class1.size <- length(C[[1]])
 class2.size <- length(C[[2]])
 axis(3, at = c(floor(class1.size/2), class1.size + floor(class2.size/2)), labels = col.classes,
  tick = FALSE, las = 1, cex.axis = 1.25, font.axis = 2, line = -1)

 return()
}
