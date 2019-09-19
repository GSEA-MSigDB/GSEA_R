#' Computes random permutation enrichment scores
#'
#' `GSEA.EnrichmentScore2` computes the weighted GSEA score of random permutations of a gene.set in gene.list
#'
#' Internal `GSEA` function.
#' Computes the weighted GSEA score of gene.set in gene.list. It is the same
#' calculation as in GSEA.EnrichmentScore but faster (x8) without producing the
#' RES, arg.RES and tag.indicator outputs.  This call is intended to be used to
#' asses the enrichment of random permutations rather than the observed one.  The
#' weighted score type is the exponent of the correlation weight: 0 (unweighted =
#' Kolmogorov-Smirnov), 1 (weighted), and 2 (over-weighted). When the score type
#' is 1 or 2 it is necessary to input the correlation vector with the values in
#' the same order as in the gene list.  Inputs: gene.list: The ordered gene list
#' (e.g. integers indicating the original position in the input dataset) gene.set:
#' A gene set (e.g. integers indicating the location of those genes in the input
#' dataset) weighted.score.type: Type of score: weight: 0 (unweighted =
#' Kolmogorov-Smirnov), 1 (weighted), and 2 (over-weighted) correl.vector: A
#' vector with the coorelations (e.g. signal to noise scores) corresponding to the
#' genes in the gene list Outputs: ES: Enrichment score (real number between -1
#' and +1)
#'
#' @keywords internal
#'

GSEA.EnrichmentScore2 <-
function(gene.list, gene.set, weighted.score.type = 1, correl.vector = NULL) {

 N <- length(gene.list)
 Nh <- length(gene.set)
 Nm <- N - Nh

 loc.vector <- vector(length = N, mode = "numeric")
 peak.res.vector <- vector(length = Nh, mode = "numeric")
 valley.res.vector <- vector(length = Nh, mode = "numeric")
 tag.correl.vector <- vector(length = Nh, mode = "numeric")
 tag.diff.vector <- vector(length = Nh, mode = "numeric")
 tag.loc.vector <- vector(length = Nh, mode = "numeric")

 loc.vector[gene.list] <- seq(1, N)
 tag.loc.vector <- loc.vector[gene.set]

 tag.loc.vector <- sort(tag.loc.vector, decreasing = F)

 if (weighted.score.type == 0) {
  tag.correl.vector <- rep(1, Nh)
 } else if (weighted.score.type == 1) {
  tag.correl.vector <- correl.vector[tag.loc.vector]
  tag.correl.vector <- abs(tag.correl.vector)
 } else if (weighted.score.type == 2) {
  tag.correl.vector <- correl.vector[tag.loc.vector] * correl.vector[tag.loc.vector]
  tag.correl.vector <- abs(tag.correl.vector)
 } else {
  tag.correl.vector <- correl.vector[tag.loc.vector]^weighted.score.type
  tag.correl.vector <- abs(tag.correl.vector)
 }

 norm.tag <- 1/sum(tag.correl.vector)
 tag.correl.vector <- tag.correl.vector * norm.tag
 norm.no.tag <- 1/Nm
 tag.diff.vector[1] <- (tag.loc.vector[1] - 1)
 tag.diff.vector[2:Nh] <- tag.loc.vector[2:Nh] - tag.loc.vector[1:(Nh - 1)] -
  1
 tag.diff.vector <- tag.diff.vector * norm.no.tag
 peak.res.vector <- cumsum(tag.correl.vector - tag.diff.vector)
 valley.res.vector <- peak.res.vector - tag.correl.vector
 max.ES <- max(peak.res.vector)
 min.ES <- min(valley.res.vector)
 ES <- signif(ifelse(max.ES > -min.ES, max.ES, min.ES), digits = 5)

 return(list(ES = ES))

}
