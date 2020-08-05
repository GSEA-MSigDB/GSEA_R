#' Run Gene Set enrichment Analysis
#'
#' `GSEA` is the main function to perform gene set enrichment analysis of a dataset
#'
#' This is a methodology for the analysis of global molecular profiles called Gene
#' Set Enrichment Analysis (GSEA). It determines whether an a priori defined set
#' of genes shows statistically significant, concordant differences between two
#' biological states (e.g. phenotypes). GSEA operates on all genes from an
#' experiment, rank ordered by the signal to noise ratio and determines whether
#' members of an a priori defined gene set are nonrandomly distributed towards the
#' top or bottom of the list and thus may correspond to an important biological
#' process. To assess significance the program uses an empirical permutation
#' procedure to test deviation from random that preserves correlations between
#' genes.  For details see Subramanian et al 2005.
#'
#' @param input.ds Input gene expression dataset file in GCT format or RNK format if preranked is specified to gsea.type
#' @param input.cls Input class vector (phenotype) file in CLS format
#' @param input.chip If collapse.dataset = TRUE, read in a CHIP formatted gene mapping file to convert dataset to gene symbols
#' @param gene.ann Depreciated parameter. Gene microarray annotation file (Affymetrix Netaffyx *.csv format) (default: none)
#' @param gs.db Gene set database in GMT format
#' @param gs.ann Depreciated parameter. Gene Set database annotation file (default: none)
#' @param output.directory Directory where to store output and results (default: .)
#' @param doc.string Documentation string used as a prefix to name result files (default: 'gsea_result')
#' @param reshuffling.type Type of permutation reshuffling: 'sample.labels' or 'gene.labels' (default: 'sample.labels')
#' @param nperm Number of random permutations (default: 1000)
#' @param weighted.score.type Enrichment correlation-based weighting: 0=no weight (KS), 1=standard weigth, 2 = over-weigth (default: 1)
#' @param nom.p.val.threshold Significance threshold for nominal p-vals for gene sets (default: -1, no thres)
#' @param fwer.p.val.threshold Significance threshold for FWER p-vals for gene sets (default: -1, no thres)
#' @param fdr.q.val.threshold Significance threshold for FDR q-vals for gene sets (default: 0.25)
#' @param topgs Besides those passing test, number of top scoring gene sets used for detailed reports (default: 20)
#' @param adjust.FDR.q.val Adjust the FDR q-vals (default: F)
#' @param gs.size.threshold.min Minimum size (in genes) for database gene sets to be considered (default: 15)
#' @param gs.size.threshold.max Maximum size (in genes) for database gene sets to be considered (default: 500)
#' @param reverse.sign Reverse direction of gene list (pos. enrichment becomes negative, etc.) (default: F)
#' @param preproc.type Preprocessing normalization: 0=none, 1=col(z-score)., 2=col(rank) and row(z-score).,
#' 3=col(rank). (default: 0)
#' @param random.seed Random number generator seed. (default: use system time)
#' @param perm.type Permutation type: 0 = unbalanced, 1 = balanced. For experts only (default: 0)
#' @param fraction Subsampling fraction. Set to 1.0 (no resampling). For experts only (default: 1.0)
#' @param replace Resampling mode (replacement or not replacement). For experts only (default: F)
#' @param collapse.dataset collapse dataset from user specified identifiers to Gene Symbols (default: FALSE)
#' @param collapse.mode Method for collapsing the dataset, accepts 'max', 'median', 'mean', 'sum', (default: NOCOLLAPSE)
#' @param save.intermediate.results save intermediate results files including ranks and permutations
#' @param use.fast.enrichment.routine If true it uses a faster GSEA.EnrichmentScore2 to compute random perm. enrichment
#' @param gsea.type Mode to run GSEA. Specify either 'GSEA' for standard mode, or 'preranked' to allow parsing of .RNK file
#' @param rank.metric Method for ranking genes. Accepts either signal-to-noise ratio 'S2N' or 'ttest' (default: S2N)
#' @return The results of the method are stored in the
#' 'output.directory' specified by the user as part of the input parameters.  The
#' results files are: - Two tab-separated global result text files (one for each
#' phenotype). These files are labeled according to the doc string prefix and the
#' phenotype name from the CLS file:
#' <doc.string>.SUMMARY.RESULTS.REPORT.<phenotype>.txt - One set of global plots.
#' They include a.- gene list correlation profile, b.- global observed and null
#' densities, c.- heat map for the entire sorted dataset, and d.- p-values vs. NES
#' plot. These plots are in a single JPEG file named
#' <doc.string>.global.plots.<phenotype>.jpg. When the program is run
#' interactively these plots appear on a window in the R GUI.  - A variable number
#' of tab-separated gene result text files according to how many sets pass any of
#' the significance thresholds ('nom.p.val.threshold,' 'fwer.p.val.threshold,' and
#' 'fdr.q.val.threshold') and how many are specified in the 'topgs' parameter.
#' These files are named: <doc.string>.<gene set name>.report.txt.  - A variable
#' number of gene set plots (one for each gene set report file). These plots
#' include a.- Gene set running enrichment 'mountain' plot, b.- gene set null
#' distribution and c.- heat map for genes in the gene set. These plots are stored
#' in a single JPEG file named <doc.string>.<gene set name>.jpg.  The format
#' (columns) for the global result files is as follows.  GS : Gene set name.  SIZE
#' : Size of the set in genes.  SOURCE : Set definition or source.  ES :
#' Enrichment score.  NES : Normalized (multiplicative rescaling) normalized
#' enrichment score.  NOM p-val : Nominal p-value (from the null distribution of
#' the gene set).  FDR q-val: False discovery rate q-values FWER p-val: Family
#' wise error rate p-values.  Tag %: Percent of gene set before running enrichment
#' peak.  Gene %: Percent of gene list before running enrichment peak.  Signal :
#' enrichment signal strength.  FDR (median): FDR q-values from the median of the
#' null distributions.  glob.p.val: P-value using a global statistic (number of
#' sets above the set's NES).  The rows are sorted by the NES values (from maximum
#' positive or negative NES to minimum) The format (columns) for the gene set
#' result files is as follows.  #: Gene number in the (sorted) gene set GENE :
#' gene name. For example the probe accession number, gene symbol or the gene
#' identifier gin the dataset.  SYMBOL : gene symbol from the gene annotation
#' file.  DESC : gene description (title) from the gene annotation file.  LIST LOC
#' : location of the gene in the sorted gene list.  S2N : signal to noise ratio
#' (correlation) of the gene in the gene list.  RES : value of the running
#' enrichment score at the gene location.  CORE_ENRICHMENT: is this gene is the
#' 'core enrichment' section of the list? Yes or No variable specifying in the
#' gene location is before (positive ES) or after (negative ES) the running
#' enrichment peak.  The rows are sorted by the gene location in the gene list.
#' The function call to GSEA returns a two element list containing the two global
#' result reports as data frames ($report1, $report2).  results1: Global output
#' report for first phenotype result2: Global putput report for second phenotype
#'
#' @examples
#' GSEA(input.ds = system.file('extdata', 'Leukemia_hgu95av2.gct', package = 'GSEA', mustWork = TRUE),
#' input.cls = system.file('extdata', 'Leukemia.cls', package = 'GSEA', mustWork = TRUE),
#' input.chip = system.file('extdata', 'Human_AFFY_HG_U95_MSigDB_7_0_final.chip',
#' package = 'GSEA', mustWork = TRUE), gs.db = system.file('extdata',
#' 'h.all.v7.0.symbols.gmt', package = 'GSEA', mustWork = TRUE),
#' collapse.dataset = TRUE, collapse.mode = 'max')
#'
#' @importFrom grDevices colors dev.cur dev.off pdf rainbow savePlot
#' @importFrom graphics axis image layout legend lines par plot points text
#' @importFrom stats density dist hclust median pnorm sd
#' @importFrom utils read.delim read.table write.table
#' @importFrom rlang .data
#' @import dplyr
#'
#' @export

GSEA <- function(input.ds, input.cls, input.chip = "NOCHIP", gene.ann = "", gs.db, 
 gs.ann = "", output.directory = getwd(), doc.string = "gsea_result", reshuffling.type = "sample.labels", 
 nperm = 1000, weighted.score.type = 1, nom.p.val.threshold = -1, fwer.p.val.threshold = -1, 
 fdr.q.val.threshold = 0.25, topgs = 20, adjust.FDR.q.val = F, gs.size.threshold.min = 15, 
 gs.size.threshold.max = 500, reverse.sign = F, preproc.type = 0, random.seed = as.integer(Sys.time()), 
 perm.type = 0, fraction = 1, replace = F, collapse.dataset = FALSE, collapse.mode = "NOCOLLAPSE", 
 save.intermediate.results = F, use.fast.enrichment.routine = T, gsea.type = "GSEA", 
 rank.metric = "S2N") {
 
 print(" *** Running Gene Set Enrichment Analysis...")
 
 # Copy input parameters to log file
 
 if (output.directory != "") {
  
  
  output.directory <- paste0(output.directory, .Platform$file.sep)
  
  filename <- paste(output.directory, doc.string, "_params.txt", sep = "", 
   collapse = "")
  
  time.string <- as.character(random.seed)
  write(paste("Run of GSEA on ", time.string), file = filename)
  
  if (is.data.frame(input.ds)) {
   # write(paste('input.ds=', quote(input.ds), sep=' '), file=filename, append=T)
  } else {
   write(paste("input.ds =", input.ds, sep = " "), file = filename, append = T)
  }
  if (gsea.type == "GSEA") {
   if (is.list(input.cls)) {
    # write(paste('input.cls=', input.cls, sep=' '), file=filename, append=T)
   } else {
    write(paste("input.cls =", input.cls, sep = " "), file = filename, 
      append = T)
   }
  }
  if (is.data.frame(gene.ann)) {
   # write(paste('gene.ann =', gene.ann, sep=' '), file=filename, append=T)
  } else {
   write(paste("gene.ann =", gene.ann, sep = " "), file = filename, append = T)
  }
  if (regexpr(pattern = ".gmt", gs.db[1]) == -1) {
   # write(paste('gs.db=', gs.db, sep=' '), file=filename, append=T)
  } else {
   write(paste("gs.db =", gs.db, sep = " "), file = filename, append = T)
  }
  if (is.data.frame(gs.ann)) {
   # write(paste('gene.ann =', gene.ann, sep=' '), file=filename, append=T)
  } else {
   write(paste("gs.ann =", gs.ann, sep = " "), file = filename, append = T)
  }
  if (gsea.type == "GSEA" & collapse.dataset == TRUE) {
   if (is.data.frame(input.chip)) {
    # write(paste('input.chip=', quote(input.chip), sep=' '), file=filename,
    # append=T)
   } else {
    write(paste("input.chip =", input.chip, sep = " "), file = filename, 
      append = T)
   }
   write(paste("collapse.mode =", collapse.mode, sep = " "), file = filename, 
    append = T)
  }
  write(paste("output.directory =", output.directory, sep = " "), file = filename, 
   append = T)
  write(paste("doc.string = ", doc.string, sep = " "), file = filename, append = T)
  write(paste("reshuffling.type =", reshuffling.type, sep = " "), file = filename, 
   append = T)
  write(paste("nperm =", nperm, sep = " "), file = filename, append = T)
  write(paste("weighted.score.type =", weighted.score.type, sep = " "), file = filename, 
   append = T)
  write(paste("nom.p.val.threshold =", nom.p.val.threshold, sep = " "), file = filename, 
   append = T)
  write(paste("fwer.p.val.threshold =", fwer.p.val.threshold, sep = " "), file = filename, 
   append = T)
  write(paste("fdr.q.val.threshold =", fdr.q.val.threshold, sep = " "), file = filename, 
   append = T)
  write(paste("topgs =", topgs, sep = " "), file = filename, append = T)
  write(paste("adjust.FDR.q.val =", adjust.FDR.q.val, sep = " "), file = filename, 
   append = T)
  write(paste("gs.size.threshold.min =", gs.size.threshold.min, sep = " "), 
   file = filename, append = T)
  write(paste("gs.size.threshold.max =", gs.size.threshold.max, sep = " "), 
   file = filename, append = T)
  write(paste("reverse.sign =", reverse.sign, sep = " "), file = filename, 
   append = T)
  write(paste("preproc.type =", preproc.type, sep = " "), file = filename, 
   append = T)
  write(paste("random.seed =", random.seed, sep = " "), file = filename, append = T)
  write(paste("perm.type =", perm.type, sep = " "), file = filename, append = T)
  if (gsea.type == "GSEA") {
   write(paste("rank.metric = ", rank.metric, sep = " "), file = filename, 
    append = T)
  } else if (gsea.type == "preranked") {
   write(paste("rank.metric = ", "preranked", sep = " "), file = filename, 
    append = T)
  }
  write(paste("fraction =", fraction, sep = " "), file = filename, append = T)
  write(paste("replace =", replace, sep = " "), file = filename, append = T)
 }
 
 # Start of GSEA methodology
 
 # Read input data matrix
 
 set.seed(seed = random.seed, kind = NULL)
 adjust.param <- 0.5
 
 gc()
 
 time1 <- proc.time()
 if (gsea.type == "GSEA") {
  if (collapse.dataset == FALSE) {
   if (is.data.frame(input.ds)) {
    dataset <- input.ds
   } else {
    dataset <- GSEA.Gct2Frame(filename = input.ds)
   }
   colnames(dataset)[1] <- "NAME"
   colnames(dataset)[2] <- "Description"
   dataset <- dataset[match(unique(dataset$NAME), dataset$NAME), ]
   dataset.ann <- dataset[, c("NAME", "Description")]
   colnames(dataset.ann)[1] <- "Gene.Symbol"
   colnames(dataset.ann)[2] <- "Gene.Title"
   rownames(dataset) <- dataset[, 1]
   gene.map <- dataset[, c(1, 2)]
   dataset <- dataset[, -1]
   dataset <- dataset[, -1]
   
  } else if (collapse.dataset == TRUE) {
   chip <- GSEA.ReadCHIPFile(file = input.chip)
   if (is.data.frame(input.ds)) {
    dataset <- input.ds
   } else {
    dataset <- GSEA.Gct2Frame(filename = input.ds)
   }
   collapseddataset <- GSEA.CollapseDataset(dataplatform = chip, gct = dataset, 
    collapse.mode = collapse.mode)
   rownames(collapseddataset) <- collapseddataset[, 1]
   gene.map <- collapseddataset[, c(1, 2)]
   collapseddataset <- collapseddataset[, -1]
   collapseddataset <- collapseddataset[, -1]
   dataset <- collapseddataset
  }

  if (rank.metric == "seq") {
   print(c("Perofrming Low Count Filtering (Preprocessing Dataset for DESeq2)"))
   dataset2 <- subset(dataset, rowSums(dataset[]) >= 10)
  }
  
  gene.labels <- row.names(dataset)
  sample.names <- colnames(dataset[1:length(colnames(dataset))])
  A <- data.matrix(dataset)
  dim(A)
  cols <- length(A[1, ])
  rows <- length(A[, 1])
  
  # preproc.type control the type of pre-processing: threshold, variation filter,
  # normalization
  
  if (preproc.type == 1) {
   # Column normalize (Z-score)
   A <- GSEA.NormalizeCols(A)
  } else if (preproc.type == 2) {
   # Column (rank) and row (Z-score) normalize column rank normalization
   for (j in 1:cols) {
    A[, j] <- rank(A[, j])
   }
   A <- GSEA.NormalizeRows(A)
  } else if (preproc.type == 3) {
   # Column (rank) norm.  column rank normalization
   for (j in 1:cols) {
    A[, j] <- rank(A[, j])
   }
  }
  
  # Read input class vector
  if (is.list(input.cls)) {
   CLS <- input.cls
  } else {
   CLS <- GSEA.ReadClsFile(file = input.cls)
  }
  class.labels <- CLS$class.v
  class.phen <- CLS$phen
  
  if (reverse.sign == T) {
   phen1 <- class.phen[2]
   phen2 <- class.phen[1]
  } else {
   phen1 <- class.phen[1]
   phen2 <- class.phen[2]
  }
  
  # sort samples according to phenotype
  
  col.index <- order(class.labels, decreasing = F)
  class.labels <- class.labels[col.index]
  sample.names <- sample.names[col.index]
  for (j in 1:rows) {
   A[j, ] <- A[j, col.index]
  }
  names(A) <- sample.names
	colnames(A) <- sample.names
 } else if (gsea.type == "preranked") {
  dataset <- read.table(input.ds, sep = "\t", header = FALSE, quote = "", stringsAsFactors = FALSE, 
   fill = TRUE)
  colnames(dataset)[1] <- "NAME"
  dataset <- dataset[match(unique(dataset$NAME), dataset$NAME), ]
  dataset.ann <- as.data.frame(dataset[, 1], stringsAsFactors = FALSE, header = FALSE)
  rownames(dataset) <- dataset[, 1]
  dataset = subset(dataset, select = -c(1))
  gene.labels <- row.names(dataset)
  dataset.ann <- cbind(dataset.ann, "NA", stringsAsFactors = FALSE)
  colnames(dataset.ann)[1] <- "Gene.Symbol"
  colnames(dataset.ann)[2] <- "Gene.Title"
  A <- data.matrix(dataset)
  cols <- length(A[1, ])
  rows <- length(A[, 1])
  phen1 <- "NA_pos"
  phen2 <- "NA_neg"
 }
 
 # Read input gene set database
 
 if (regexpr(pattern = ".gmt", gs.db[1]) == -1) {
  temp <- gs.db
 } else {
  temp <- readLines(gs.db)
 }
 
 max.Ng <- length(temp)
 temp.size.G <- vector(length = max.Ng, mode = "numeric")
 for (i in 1:max.Ng) {
  temp.size.G[i] <- length(unlist(strsplit(temp[[i]], "\\t"))) - 2
 }
 
 max.size.G <- max(temp.size.G)
 min.size.G <- min(temp.size.G)
 gs <- matrix(rep("null", max.Ng * max.size.G), nrow = max.Ng, ncol = max.size.G)
 temp.names <- vector(length = max.Ng, mode = "character")
 temp.desc <- vector(length = max.Ng, mode = "character")
 gs.count <- 1
 for (i in 1:max.Ng) {
  gene.set.size <- length(unlist(strsplit(temp[[i]], "\\t"))) - 2
  gs.line <- noquote(unlist(strsplit(temp[[i]], "\\t")))
  gene.set.name <- gs.line[1]
  gene.set.desc <- gs.line[2]
  gene.set.tags <- vector(length = gene.set.size, mode = "character")
  for (j in 1:gene.set.size) {
   gene.set.tags[j] <- gs.line[j + 2]
  }
  existing.set <- is.element(gene.set.tags, gene.labels)
  set.size <- length(existing.set[existing.set == T])
  if ((set.size < gs.size.threshold.min) || (set.size > gs.size.threshold.max)) 
   next
  temp.size.G[gs.count] <- set.size
  gs[gs.count, ] <- c(gene.set.tags[existing.set], rep("null", max.size.G - 
   temp.size.G[gs.count]))
  temp.names[gs.count] <- gene.set.name
  temp.desc[gs.count] <- gene.set.desc
  gs.count <- gs.count + 1
 }
 Ng <- gs.count - 1
 gs.names <- vector(length = Ng, mode = "character")
 gs.desc <- vector(length = Ng, mode = "character")
 size.G <- vector(length = Ng, mode = "numeric")
 gs.names <- temp.names[1:Ng]
 gs.desc <- temp.desc[1:Ng]
 size.G <- temp.size.G[1:Ng]
 
 N <- length(A[, 1])
 Ns <- length(A[1, ])
 
 print(c("Number of genes:", N))
 print(c("Number of Gene Sets:", Ng))
 print(c("Number of samples:", Ns))
 print(c("Original number of Gene Sets:", max.Ng))
 print(c("Maximum gene set size:", max.size.G))
 print(c("Minimum gene set size:", min.size.G))
 
 # Read gene and gene set annotations if gene annotation file was provided
 
 all.gene.descs <- vector(length = N, mode = "character")
 all.gene.symbols <- vector(length = N, mode = "character")
 all.gs.descs <- vector(length = Ng, mode = "character")
 
 if (exists("chip") == TRUE) {
  gene.ann <- unique(chip[, c("Gene.Symbol", "Gene.Title")])
  print(c("Replacing gene descriptions from input dataset with annotations from mapping CHIP."))
 }
 if (exists("chip") == FALSE) {
  gene.ann <- dataset.ann
 }
 
 if (is.data.frame(gene.ann)) {
  temp <- gene.ann
  a.size <- length(temp[, 1])
  print(paste("Number of gene annotation file entries:", a.size))
  accs <- as.character(temp[, 1])
  locs <- match(gene.labels, accs)
  all.gene.descs <- as.character(temp[locs, "Gene.Title"])
  all.gene.symbols <- as.character(temp[locs, "Gene.Symbol"])
  rm(temp)
 } else if (gene.ann == "") {
  
  for (i in 1:N) {
   all.gene.descs[i] <- gene.map[i, 2]
   all.gene.symbols[i] <- gene.labels[i]
  }
  
 }
 
 if (is.data.frame(gs.ann)) {
  temp <- gs.ann
  a.size <- length(temp[, 1])
  print(paste("Number of gene set annotation file entries:", a.size))
  accs <- as.character(temp[, 1])
  locs <- match(gs.names, accs)
  all.gs.descs <- as.character(temp[locs, "SOURCE"])
  rm(temp)
 } else if (gs.ann == "") {
  for (i in 1:Ng) {
   all.gs.descs[i] <- gs.desc[i]
  }
 } else {
  temp <- read.delim(gs.ann, header = T, sep = "\t", comment.char = "", as.is = T)
  a.size <- length(temp[, 1])
  print(paste("Number of gene set annotation file entries:", a.size))
  accs <- as.character(temp[, 1])
  locs <- match(gs.names, accs)
  all.gs.descs <- as.character(temp[locs, "SOURCE"])
  rm(temp)
 }
 
 
 Obs.indicator <- matrix(nrow = Ng, ncol = N)
 Obs.RES <- matrix(nrow = Ng, ncol = N)
 
 Obs.ES <- vector(length = Ng, mode = "numeric")
 Obs.arg.ES <- vector(length = Ng, mode = "numeric")
 Obs.ES.norm <- vector(length = Ng, mode = "numeric")
 
 time2 <- proc.time()
 
 # GSEA methodology
 
 # Compute observed and random permutation gene rankings
 
 obs.rnk <- vector(length = N, mode = "numeric")
 signal.strength <- vector(length = Ng, mode = "numeric")
 tag.frac <- vector(length = Ng, mode = "numeric")
 gene.frac <- vector(length = Ng, mode = "numeric")
 coherence.ratio <- vector(length = Ng, mode = "numeric")
 obs.phi.norm <- matrix(nrow = Ng, ncol = nperm)
 correl.matrix <- matrix(nrow = N, ncol = nperm)
 obs.correl.matrix <- matrix(nrow = N, ncol = nperm)
 order.matrix <- matrix(nrow = N, ncol = nperm)
 obs.order.matrix <- matrix(nrow = N, ncol = nperm)
 
 nperm.per.call <- 100
 n.groups <- nperm%/%nperm.per.call
 n.rem <- nperm%%nperm.per.call
 n.perms <- c(rep(nperm.per.call, n.groups), n.rem)
 n.ends <- cumsum(n.perms)
 n.starts <- n.ends - n.perms + 1
 
 if (n.rem == 0) {
  n.tot <- n.groups
 } else {
  n.tot <- n.groups + 1
 }
 if (gsea.type == "GSEA") {
  for (nk in 1:n.tot) {
   call.nperm <- n.perms[nk]
   
   print(paste("Computing ranked list for actual and permuted phenotypes.......permutations: ", 
    n.starts[nk], "--", n.ends[nk], sep = " "))
   
   O <- GSEA.GeneRanking(A, class.labels, gene.labels, call.nperm, permutation.type = perm.type, 
    sigma.correction = "GeneCluster", fraction = fraction, replace = replace, 
    reverse.sign = reverse.sign, rank.metric)
   gc()
   if (length(O$obs.rnk.matrix[is.na(O$obs.rnk.matrix)])>0) {
    warning(print(length(O$obs.rnk.matrix[is.na(O$obs.rnk.matrix)]))," N/A values were found in the ranked list. Setting N/As to Zero because thse cause GSEA to fail.")
    O$obs.rnk.matrix[is.na(O$obs.rnk.matrix)]<-0
   }
   if (length(O$rnk.matrix[is.na(O$rnk.matrix)])>0) {
    warning(print(length(O$rnk.matrix[is.na(O$rnk.matrix)])), " N/A values were found in the ranked list. Setting N/As to Zero because thse cause GSEA to fail.")
    O$rnk.matrix[is.na(O$rnk.matrix)]<-0
   }
   order.matrix[, n.starts[nk]:n.ends[nk]] <- O$order.matrix
   obs.order.matrix[, n.starts[nk]:n.ends[nk]] <- O$obs.order.matrix
   correl.matrix[, n.starts[nk]:n.ends[nk]] <- O$rnk.matrix
   obs.correl.matrix[, n.starts[nk]:n.ends[nk]] <- O$obs.rnk.matrix
   rm(O)
  }
  obs.rnk <- apply(obs.correl.matrix, 1, median)  # using median to assign enrichment scores
  obs.index <- order(obs.rnk, decreasing = T)
  obs.rnk <- sort(obs.rnk, decreasing = T)
  
  obs.gene.labels <- gene.labels[obs.index]
  obs.gene.descs <- all.gene.descs[obs.index]
  obs.gene.symbols <- all.gene.symbols[obs.index]
 } else if (gsea.type == "preranked") {
  print(paste("Skipping calculating gene rankings... using pre-ranked list."))
  obs.rnk <- unname(A[, 1])
  obs.index <- order(obs.rnk, decreasing = T)
  obs.rnk <- sort(obs.rnk, decreasing = T)
  gene.list2 <- obs.index
 }
 
 for (r in 1:nperm) {
  correl.matrix[, r] <- correl.matrix[order.matrix[, r], r]
 }
 for (r in 1:nperm) {
  obs.correl.matrix[, r] <- obs.correl.matrix[obs.order.matrix[, r], r]
 }
 if (gsea.type == "preranked") {
  obs.gene.labels <- gene.labels[obs.index]
  obs.gene.descs <- obs.gene.labels
  obs.gene.symbols <- obs.gene.labels
  obs.order.matrix[, 1:1000] <- gene.list2
  obs.correl.matrix[, 1:1000] <- obs.rnk
 }
 
 gene.list2 <- obs.index
 for (i in 1:Ng) {
 try{
  print(paste("Computing observed enrichment for gene set:", i, gs.names[i], 
   sep = " "))
  gene.set <- gs[i, gs[i, ] != "null"]
  gene.set2 <- vector(length = length(gene.set), mode = "numeric")
  gene.set2 <- match(gene.set, gene.labels)
  GSEA.results <- GSEA.EnrichmentScore(gene.list = gene.list2, gene.set = gene.set2, 
   weighted.score.type = weighted.score.type, correl.vector = obs.rnk)
  Obs.ES[i] <- GSEA.results$ES
  Obs.arg.ES[i] <- GSEA.results$arg.ES
  Obs.RES[i, ] <- GSEA.results$RES
  Obs.indicator[i, ] <- GSEA.results$indicator
  if (Obs.ES[i] >= 0) {
   # compute signal strength
   tag.frac[i] <- sum(Obs.indicator[i, 1:Obs.arg.ES[i]])/size.G[i]
   gene.frac[i] <- Obs.arg.ES[i]/N
  } else {
   tag.frac[i] <- sum(Obs.indicator[i, Obs.arg.ES[i]:N])/size.G[i]
   gene.frac[i] <- (N - Obs.arg.ES[i] + 1)/N
  }
  signal.strength[i] <- tag.frac[i] * (1 - gene.frac[i]) * (N/(N - size.G[i]))
 }
 }
 
 # Compute enrichment for random permutations
 
 phi <- matrix(nrow = Ng, ncol = nperm)
 phi.norm <- matrix(nrow = Ng, ncol = nperm)
 obs.phi <- matrix(nrow = Ng, ncol = nperm)
 
 if (reshuffling.type == "sample.labels") {
  # reshuffling phenotype labels
  for (i in 1:Ng) {
   try {
   print(paste("Computing random permutations' enrichment for gene set:", 
    i, gs.names[i], sep = " "))
   gene.set <- gs[i, gs[i, ] != "null"]
   gene.set2 <- vector(length = length(gene.set), mode = "numeric")
   gene.set2 <- match(gene.set, gene.labels)
   for (r in 1:nperm) {
    gene.list2 <- order.matrix[, r]
    if (use.fast.enrichment.routine == F) {
      GSEA.results <- GSEA.EnrichmentScore(gene.list = gene.list2, gene.set = gene.set2, 
     weighted.score.type = weighted.score.type, correl.vector = correl.matrix[, 
       r])
    } else {
      GSEA.results <- GSEA.EnrichmentScore2(gene.list = gene.list2, gene.set = gene.set2, 
     weighted.score.type = weighted.score.type, correl.vector = correl.matrix[, 
       r])
    }
    phi[i, r] <- GSEA.results$ES
   }
   if (fraction < 1) {
    # if resampling then compute ES for all observed rankings
    for (r in 1:nperm) {
      obs.gene.list2 <- obs.order.matrix[, r]
      if (use.fast.enrichment.routine == F) {
     GSEA.results <- GSEA.EnrichmentScore(gene.list = obs.gene.list2, 
       gene.set = gene.set2, weighted.score.type = weighted.score.type, 
       correl.vector = obs.correl.matrix[, r])
      } else {
     GSEA.results <- GSEA.EnrichmentScore2(gene.list = obs.gene.list2, 
       gene.set = gene.set2, weighted.score.type = weighted.score.type, 
       correl.vector = obs.correl.matrix[, r])
      }
      obs.phi[i, r] <- GSEA.results$ES
    }
   } else {
    # if no resampling then compute only one column (and fill the others with the
    # same value)
    obs.gene.list2 <- obs.order.matrix[, 1]
    if (use.fast.enrichment.routine == F) {
      GSEA.results <- GSEA.EnrichmentScore(gene.list = obs.gene.list2, 
     gene.set = gene.set2, weighted.score.type = weighted.score.type, 
     correl.vector = obs.correl.matrix[, r])
    } else {
      GSEA.results <- GSEA.EnrichmentScore2(gene.list = obs.gene.list2, 
     gene.set = gene.set2, weighted.score.type = weighted.score.type, 
     correl.vector = obs.correl.matrix[, r])
    }
    obs.phi[i, 1] <- GSEA.results$ES
    for (r in 2:nperm) {
      obs.phi[i, r] <- obs.phi[i, 1]
    }
   }
   gc()
   }
  }
  
 } else if (reshuffling.type == "gene.labels") {
  # reshuffling gene labels
  for (i in 1:Ng) {
   try{
   print(paste("Computing random permutations' enrichment for gene set:", 
    i, gs.names[i], sep = " "))
   gene.set <- gs[i, gs[i, ] != "null"]
   gene.set2 <- vector(length = length(gene.set), mode = "numeric")
   gene.set2 <- match(gene.set, gene.labels)
   for (r in 1:nperm) {
    reshuffled.gene.labels <- sample(1:rows)
    if (use.fast.enrichment.routine == F) {
      GSEA.results <- GSEA.EnrichmentScore(gene.list = reshuffled.gene.labels, 
     gene.set = gene.set2, weighted.score.type = weighted.score.type, 
     correl.vector = obs.rnk)
    } else {
      GSEA.results <- GSEA.EnrichmentScore2(gene.list = reshuffled.gene.labels, 
     gene.set = gene.set2, weighted.score.type = weighted.score.type, 
     correl.vector = obs.rnk)
    }
    phi[i, r] <- GSEA.results$ES
   }
   if (fraction < 1) {
    # if resampling then compute ES for all observed rankings
    for (r in 1:nperm) {
      obs.gene.list2 <- obs.order.matrix[, r]
      if (use.fast.enrichment.routine == F) {
     GSEA.results <- GSEA.EnrichmentScore(gene.list = obs.gene.list2, 
       gene.set = gene.set2, weighted.score.type = weighted.score.type, 
       correl.vector = obs.correl.matrix[, r])
      } else {
     GSEA.results <- GSEA.EnrichmentScore2(gene.list = obs.gene.list2, 
       gene.set = gene.set2, weighted.score.type = weighted.score.type, 
       correl.vector = obs.correl.matrix[, r])
      }
      obs.phi[i, r] <- GSEA.results$ES
    }
   } else {
    # if no resampling then compute only one column (and fill the others with the
    # same value)
    obs.gene.list2 <- obs.order.matrix[, 1]
    if (use.fast.enrichment.routine == F) {
      GSEA.results <- GSEA.EnrichmentScore(gene.list = obs.gene.list2, 
     gene.set = gene.set2, weighted.score.type = weighted.score.type, 
     correl.vector = obs.correl.matrix[, r])
    } else {
      GSEA.results <- GSEA.EnrichmentScore2(gene.list = obs.gene.list2, 
     gene.set = gene.set2, weighted.score.type = weighted.score.type, 
     correl.vector = obs.correl.matrix[, r])
    }
    obs.phi[i, 1] <- GSEA.results$ES
    for (r in 2:nperm) {
      obs.phi[i, r] <- obs.phi[i, 1]
    }
   }
   gc()
   }
  }
 }
 
 # Compute 3 types of p-values
 
 # Find nominal p-values
 
 print("Computing nominal p-values...")
 
 p.vals <- matrix(0, nrow = Ng, ncol = 2)
 
 for (i in 1:Ng) {
  pos.phi <- NULL
  neg.phi <- NULL
  for (j in 1:nperm) {
   if (phi[i, j] >= 0) {
    pos.phi <- c(pos.phi, phi[i, j])
   } else {
    neg.phi <- c(neg.phi, phi[i, j])
   }
  }
  ES.value <- Obs.ES[i]
  if (ES.value >= 0) {
   p.vals[i, 1] <- signif(sum(pos.phi >= ES.value)/length(pos.phi), digits = 5)
  } else {
   p.vals[i, 1] <- signif(sum(neg.phi <= ES.value)/length(neg.phi), digits = 5)
  }
 }
 
 # Find effective size
 
 erf <- function(x) {
  2 * pnorm(sqrt(2) * x)
 }
 
 KS.mean <- function(N) {
  # KS mean as a function of set size N
  S <- 0
  for (k in -100:100) {
   if (k == 0) 
    next
   S <- S + 4 * (-1)^(k + 1) * (0.25 * exp(-2 * k * k * N) - sqrt(2 * pi) * 
    erf(sqrt(2 * N) * k)/(16 * k * sqrt(N)))
  }
  return(abs(S))
 }
 
 # KS.mean.table <- vector(length=5000, mode='numeric')
 
 # for (i in 1:5000) { KS.mean.table[i] <- KS.mean(i) }
 
 # KS.size <- vector(length=Ng, mode='numeric')
 
 # Rescaling normalization for each gene set null
 
 print("Computing rescaling normalization for each gene set null...")
 
 for (i in 1:Ng) {
  pos.phi <- NULL
  neg.phi <- NULL
  for (j in 1:nperm) {
   if (phi[i, j] >= 0) {
    pos.phi <- c(pos.phi, phi[i, j])
   } else {
    neg.phi <- c(neg.phi, phi[i, j])
   }
  }
  pos.m <- mean(pos.phi)
  neg.m <- mean(abs(neg.phi))
  
  # if (Obs.ES[i] >= 0) { KS.size[i] <- which.min(abs(KS.mean.table - pos.m)) }
  # else { KS.size[i] <- which.min(abs(KS.mean.table - neg.m)) }
  
  pos.phi <- pos.phi/pos.m
  neg.phi <- neg.phi/neg.m
  for (j in 1:nperm) {
   if (phi[i, j] >= 0) {
    phi.norm[i, j] <- phi[i, j]/pos.m
   } else {
    phi.norm[i, j] <- phi[i, j]/neg.m
   }
  }
  for (j in 1:nperm) {
   if (obs.phi[i, j] >= 0) {
    obs.phi.norm[i, j] <- obs.phi[i, j]/pos.m
   } else {
    obs.phi.norm[i, j] <- obs.phi[i, j]/neg.m
   }
  }
  if (Obs.ES[i] >= 0) {
   Obs.ES.norm[i] <- Obs.ES[i]/pos.m
  } else {
   Obs.ES.norm[i] <- Obs.ES[i]/neg.m
  }
 }
 
 
 # Save intermedite results
 
 if (save.intermediate.results == T) {
  
  filename <- paste(output.directory, doc.string, ".phi.txt", sep = "", collapse = "")
  write.table(phi, file = filename, quote = F, col.names = F, row.names = F, 
   sep = "\t")
  
  filename <- paste(output.directory, doc.string, ".obs.phi.txt", sep = "", 
   collapse = "")
  write.table(obs.phi, file = filename, quote = F, col.names = F, row.names = F, 
   sep = "\t")
  
  filename <- paste(output.directory, doc.string, ".phi.norm.txt", sep = "", 
   collapse = "")
  write.table(phi.norm, file = filename, quote = F, col.names = F, row.names = F, 
   sep = "\t")
  
  filename <- paste(output.directory, doc.string, ".obs.phi.norm.txt", sep = "", 
   collapse = "")
  write.table(obs.phi.norm, file = filename, quote = F, col.names = F, row.names = F, 
   sep = "\t")
  
  filename <- paste(output.directory, doc.string, ".Obs.ES.txt", sep = "", 
   collapse = "")
  write.table(Obs.ES, file = filename, quote = F, col.names = F, row.names = F, 
   sep = "\t")
  
  filename <- paste(output.directory, doc.string, ".Obs.ES.norm.txt", sep = "", 
   collapse = "")
  write.table(Obs.ES.norm, file = filename, quote = F, col.names = F, row.names = F, 
   sep = "\t")
 }
 
 # Compute FWER p-vals
 
 print("Computing FWER p-values...")
 
 max.ES.vals.p <- NULL
 max.ES.vals.n <- NULL
 for (j in 1:nperm) {
  pos.phi <- NULL
  neg.phi <- NULL
  for (i in 1:Ng) {
   if (phi.norm[i, j] >= 0) {
    pos.phi <- c(pos.phi, phi.norm[i, j])
   } else {
    neg.phi <- c(neg.phi, phi.norm[i, j])
   }
  }
  if (length(pos.phi) > 0) {
   max.ES.vals.p <- c(max.ES.vals.p, max(pos.phi))
  }
  if (length(neg.phi) > 0) {
   max.ES.vals.n <- c(max.ES.vals.n, min(neg.phi))
  }
 }
 for (i in 1:Ng) {
  ES.value <- Obs.ES.norm[i]
  if (Obs.ES.norm[i] >= 0) {
   p.vals[i, 2] <- signif(sum(max.ES.vals.p >= ES.value)/length(max.ES.vals.p), 
    digits = 5)
  } else {
   p.vals[i, 2] <- signif(sum(max.ES.vals.n <= ES.value)/length(max.ES.vals.n), 
    digits = 5)
  }
 }
 
 
 # Compute FDRs
 
 print("Computing FDR q-values...")
 
 NES <- vector(length = Ng, mode = "numeric")
 phi.norm.mean <- vector(length = Ng, mode = "numeric")
 obs.phi.norm.mean <- vector(length = Ng, mode = "numeric")
 phi.norm.median <- vector(length = Ng, mode = "numeric")
 obs.phi.norm.median <- vector(length = Ng, mode = "numeric")
 phi.norm.mean <- vector(length = Ng, mode = "numeric")
 obs.phi.mean <- vector(length = Ng, mode = "numeric")
 FDR.mean <- vector(length = Ng, mode = "numeric")
 FDR.median <- vector(length = Ng, mode = "numeric")
 phi.norm.median.d <- vector(length = Ng, mode = "numeric")
 obs.phi.norm.median.d <- vector(length = Ng, mode = "numeric")
 
 Obs.ES.index <- order(Obs.ES.norm, decreasing = T)
 Orig.index <- seq(1, Ng)
 Orig.index <- Orig.index[Obs.ES.index]
 Orig.index <- order(Orig.index, decreasing = F)
 Obs.ES.norm.sorted <- Obs.ES.norm[Obs.ES.index]
 gs.names.sorted <- gs.names[Obs.ES.index]
 
 for (k in 1:Ng) {
  NES[k] <- Obs.ES.norm.sorted[k]
  ES.value <- NES[k]
  count.col <- vector(length = nperm, mode = "numeric")
  obs.count.col <- vector(length = nperm, mode = "numeric")
  for (i in 1:nperm) {
   phi.vec <- phi.norm[, i]
   obs.phi.vec <- obs.phi.norm[, i]
   if (ES.value >= 0) {
    count.col.norm <- sum(phi.vec >= 0)
    obs.count.col.norm <- sum(obs.phi.vec >= 0)
    count.col[i] <- ifelse(count.col.norm > 0, sum(phi.vec >= ES.value)/count.col.norm, 
      0)
    obs.count.col[i] <- ifelse(obs.count.col.norm > 0, sum(obs.phi.vec >= 
      ES.value)/obs.count.col.norm, 0)
   } else {
    count.col.norm <- sum(phi.vec < 0)
    obs.count.col.norm <- sum(obs.phi.vec < 0)
    count.col[i] <- ifelse(count.col.norm > 0, sum(phi.vec <= ES.value)/count.col.norm, 
      0)
    obs.count.col[i] <- ifelse(obs.count.col.norm > 0, sum(obs.phi.vec <= 
      ES.value)/obs.count.col.norm, 0)
   }
  }
  phi.norm.mean[k] <- mean(count.col)
  obs.phi.norm.mean[k] <- mean(obs.count.col)
  phi.norm.median[k] <- median(count.col)
  obs.phi.norm.median[k] <- median(obs.count.col)
  FDR.mean[k] <- ifelse(phi.norm.mean[k]/obs.phi.norm.mean[k] < 1, phi.norm.mean[k]/obs.phi.norm.mean[k], 
   1)
  FDR.median[k] <- ifelse(phi.norm.median[k]/obs.phi.norm.median[k] < 1, phi.norm.median[k]/obs.phi.norm.median[k], 
   1)
 }
 
 # adjust q-values
 
 if (adjust.FDR.q.val == T) {
  pos.nes <- length(NES[NES >= 0])
  min.FDR.mean <- FDR.mean[pos.nes]
  min.FDR.median <- FDR.median[pos.nes]
  for (k in seq(pos.nes - 1, 1, -1)) {
   if (FDR.mean[k] < min.FDR.mean) {
    min.FDR.mean <- FDR.mean[k]
   }
   if (min.FDR.mean < FDR.mean[k]) {
    FDR.mean[k] <- min.FDR.mean
   }
  }
  
  neg.nes <- pos.nes + 1
  min.FDR.mean <- FDR.mean[neg.nes]
  min.FDR.median <- FDR.median[neg.nes]
  for (k in seq(neg.nes + 1, Ng)) {
   if (FDR.mean[k] < min.FDR.mean) {
    min.FDR.mean <- FDR.mean[k]
   }
   if (min.FDR.mean < FDR.mean[k]) {
    FDR.mean[k] <- min.FDR.mean
   }
  }
 }
 
 obs.phi.norm.mean.sorted <- obs.phi.norm.mean[Orig.index]
 phi.norm.mean.sorted <- phi.norm.mean[Orig.index]
 FDR.mean.sorted <- FDR.mean[Orig.index]
 FDR.median.sorted <- FDR.median[Orig.index]
 
 # Compute global statistic
 
 glob.p.vals <- vector(length = Ng, mode = "numeric")
 NULL.pass <- vector(length = nperm, mode = "numeric")
 OBS.pass <- vector(length = nperm, mode = "numeric")
 
 for (k in 1:Ng) {
  NES[k] <- Obs.ES.norm.sorted[k]
  if (NES[k] >= 0) {
   for (i in 1:nperm) {
    NULL.pos <- sum(phi.norm[, i] >= 0)
    NULL.pass[i] <- ifelse(NULL.pos > 0, sum(phi.norm[, i] >= NES[k])/NULL.pos, 
      0)
    OBS.pos <- sum(obs.phi.norm[, i] >= 0)
    OBS.pass[i] <- ifelse(OBS.pos > 0, sum(obs.phi.norm[, i] >= NES[k])/OBS.pos, 
      0)
   }
  } else {
   for (i in 1:nperm) {
    NULL.neg <- sum(phi.norm[, i] < 0)
    NULL.pass[i] <- ifelse(NULL.neg > 0, sum(phi.norm[, i] <= NES[k])/NULL.neg, 
      0)
    OBS.neg <- sum(obs.phi.norm[, i] < 0)
    OBS.pass[i] <- ifelse(OBS.neg > 0, sum(obs.phi.norm[, i] <= NES[k])/OBS.neg, 
      0)
   }
  }
  glob.p.vals[k] <- sum(NULL.pass >= mean(OBS.pass))/nperm
 }
 glob.p.vals.sorted <- glob.p.vals[Orig.index]
 
 # Produce results report
 
 print("Producing result tables and plots...")
 
 Obs.ES <- signif(Obs.ES, digits = 5)
 Obs.ES.norm <- signif(Obs.ES.norm, digits = 5)
 p.vals <- signif(p.vals, digits = nchar(nperm))
 signal.strength <- signif(signal.strength, digits = 3)
 tag.frac <- signif(tag.frac, digits = 3)
 gene.frac <- signif(gene.frac, digits = 3)
 FDR.mean.sorted <- signif(FDR.mean.sorted, digits = 5)
 FDR.median.sorted <- signif(FDR.median.sorted, digits = 5)
 glob.p.vals.sorted <- signif(glob.p.vals.sorted, digits = 5)
 
 report <- data.frame(cbind(gs.names, size.G, all.gs.descs, Obs.ES, Obs.ES.norm, 
  p.vals[, 1], FDR.mean.sorted, p.vals[, 2], tag.frac, gene.frac, signal.strength, 
  FDR.median.sorted, glob.p.vals.sorted), stringsAsFactors = FALSE)
 names(report) <- c("GS", "SIZE", "SOURCE", "ES", "NES", "NOM p-val", "FDR q-val", 
  "FWER p-val", "Tag %", "Gene %", "Signal", "FDR (median)", "glob.p.val")
 # print(report)
 report2 <- report
 report.index2 <- order(Obs.ES.norm, decreasing = T)
 for (i in 1:Ng) {
  report2[i, ] <- report[report.index2[i], ]
 }
 report3 <- report
 report.index3 <- order(Obs.ES.norm, decreasing = F)
 for (i in 1:Ng) {
  report3[i, ] <- report[report.index3[i], ]
 }
 phen1.rows <- length(Obs.ES.norm[Obs.ES.norm >= 0])
 phen2.rows <- length(Obs.ES.norm[Obs.ES.norm < 0])
 report.phen1 <- report2[1:phen1.rows, ]
 report.phen2 <- report3[1:phen2.rows, ]
 
 if (output.directory != "") {
  if (phen1.rows > 0) {
   filename <- paste(output.directory, doc.string, ".SUMMARY.RESULTS.REPORT.", 
    phen1, ".txt", sep = "", collapse = "")
   write.table(report.phen1, file = filename, quote = F, row.names = F, 
    sep = "\t")
  }
  if (phen2.rows > 0) {
   filename <- paste(output.directory, doc.string, ".SUMMARY.RESULTS.REPORT.", 
    phen2, ".txt", sep = "", collapse = "")
   write.table(report.phen2, file = filename, quote = F, row.names = F, 
    sep = "\t")
  }
 }
 
 # Global plots
 
 if (output.directory != "") {
  glob.filename <- paste(output.directory, doc.string, ".global.plots.pdf", 
   sep = "", collapse = "")
  pdf(file = glob.filename, height = 10, width = 10)
 }
 
 nf <- layout(matrix(c(1, 2, 3, 4), 2, 2, byrow = T), c(1, 1), c(1, 1), TRUE)
 
 # plot S2N correlation profile
 
 location <- 1:N
 max.corr <- max(obs.rnk)
 min.corr <- min(obs.rnk)
 
 if (gsea.type == "GSEA") {
  if (rank.metric == "S2N") {
   x <- plot(location, obs.rnk, ylab = "Signal to Noise Ratio (S2N)", xlab = "Gene List Location", 
    main = "Gene List Correlation (S2N) Profile", type = "l", lwd = 2, 
    cex = 0.9, col = 1)
  }
  if (rank.metric == "ttest") {
   x <- plot(location, obs.rnk, ylab = "T-Test", xlab = "Gene List Location", 
    main = "Gene List Correlation (T-Test) Profile", type = "l", lwd = 2, 
    cex = 0.9, col = 1)
  }
  if (rank.metric == "seq") {
   x <- plot(location, obs.rnk, ylab = "DESeq2 Log2(FC)", xlab = "Gene List Location", 
    main = "Gene List Correlation (Log2(FC)) Profile", type = "l", lwd = 2, 
    cex = 0.9, col = 1)
  }
 } else if (gsea.type == "preranked") {
  x <- plot(location, obs.rnk, ylab = "User Rank Metric", xlab = "Gene List Location", 
   main = "Gene List Correlation Profile", type = "l", lwd = 2, cex = 0.9, 
   col = 1)
 }
 
 for (i in seq(1, N, 20)) {
  lines(c(i, i), c(0, obs.rnk[i]), lwd = 3, cex = 0.9, col = colors()[12])  # shading of correlation plot
 }
 x <- points(location, obs.rnk, type = "l", lwd = 2, cex = 0.9, col = 1)
 lines(c(1, N), c(0, 0), lwd = 2, lty = 1, cex = 0.9, col = 1)  # zero correlation horizontal line
 temp <- order(abs(obs.rnk), decreasing = T)
 arg.correl <- temp[N]
 lines(c(arg.correl, arg.correl), c(min.corr, 0.7 * max.corr), lwd = 2, lty = 3, 
  cex = 0.9, col = 1)  # zero correlation vertical line
 
 area.bias <- signif(100 * (sum(obs.rnk[1:arg.correl]) + sum(obs.rnk[arg.correl:N]))/sum(abs(obs.rnk[1:N])), 
  digits = 3)
 area.phen <- ifelse(area.bias >= 0, phen1, phen2)
 delta.string <- paste("Corr. Area Bias to \"", area.phen, "\" =", abs(area.bias), 
  "%", sep = "", collapse = "")
 zero.crossing.string <- paste("Zero Crossing at location ", arg.correl, " (", 
  signif(100 * arg.correl/N, digits = 3), " %)")
 leg.txt <- c(delta.string, zero.crossing.string)
 legend(x = N/10, y = max.corr, bty = "n", bg = "white", legend = leg.txt, cex = 0.9)
 
 leg.txt <- paste("\"", phen1, "\" ", sep = "", collapse = "")
 text(x = 1, y = -0.05 * max.corr, adj = c(0, 1), labels = leg.txt, cex = 0.9)
 
 leg.txt <- paste("\"", phen2, "\" ", sep = "", collapse = "")
 text(x = N, y = 0.05 * max.corr, adj = c(1, 0), labels = leg.txt, cex = 0.9)
 
 if (Ng > 1) 
  {
   # make these plots only if there are multiple gene sets.
   
   # compute plots of actual (weighted) null and observed
   
   phi.densities.pos <- matrix(0, nrow = 512, ncol = nperm)
   phi.densities.neg <- matrix(0, nrow = 512, ncol = nperm)
   obs.phi.densities.pos <- matrix(0, nrow = 512, ncol = nperm)
   obs.phi.densities.neg <- matrix(0, nrow = 512, ncol = nperm)
   phi.density.mean.pos <- vector(length = 512, mode = "numeric")
   phi.density.mean.neg <- vector(length = 512, mode = "numeric")
   obs.phi.density.mean.pos <- vector(length = 512, mode = "numeric")
   obs.phi.density.mean.neg <- vector(length = 512, mode = "numeric")
   phi.density.median.pos <- vector(length = 512, mode = "numeric")
   phi.density.median.neg <- vector(length = 512, mode = "numeric")
   obs.phi.density.median.pos <- vector(length = 512, mode = "numeric")
   obs.phi.density.median.neg <- vector(length = 512, mode = "numeric")
   x.coor.pos <- vector(length = 512, mode = "numeric")
   x.coor.neg <- vector(length = 512, mode = "numeric")
   
   for (i in 1:nperm) {
    pos.phi <- phi.norm[phi.norm[, i] >= 0, i]
    if (length(pos.phi) > 2) {
      temp <- density(pos.phi, adjust = adjust.param, n = 512, from = 0, 
     to = 3.5)
    } else {
      temp <- list(x = 3.5 * (seq(1, 512) - 1)/512, y = rep(0.001, 512))
    }
    phi.densities.pos[, i] <- temp$y
    norm.factor <- sum(phi.densities.pos[, i])
    phi.densities.pos[, i] <- phi.densities.pos[, i]/norm.factor
    if (i == 1) {
      x.coor.pos <- temp$x
    }
    
    neg.phi <- phi.norm[phi.norm[, i] < 0, i]
    if (length(neg.phi) > 2) {
      temp <- density(neg.phi, adjust = adjust.param, n = 512, from = -3.5, 
     to = 0)
    } else {
      temp <- list(x = 3.5 * (seq(1, 512) - 1)/512, y = rep(0.001, 512))
    }
    phi.densities.neg[, i] <- temp$y
    norm.factor <- sum(phi.densities.neg[, i])
    phi.densities.neg[, i] <- phi.densities.neg[, i]/norm.factor
    if (i == 1) {
      x.coor.neg <- temp$x
    }
    pos.phi <- obs.phi.norm[obs.phi.norm[, i] >= 0, i]
    if (length(pos.phi) > 2) {
      temp <- density(pos.phi, adjust = adjust.param, n = 512, from = 0, 
     to = 3.5)
    } else {
      temp <- list(x = 3.5 * (seq(1, 512) - 1)/512, y = rep(0.001, 512))
    }
    obs.phi.densities.pos[, i] <- temp$y
    norm.factor <- sum(obs.phi.densities.pos[, i])
    obs.phi.densities.pos[, i] <- obs.phi.densities.pos[, i]/norm.factor
    
    neg.phi <- obs.phi.norm[obs.phi.norm[, i] < 0, i]
    if (length(neg.phi) > 2) {
      temp <- density(neg.phi, adjust = adjust.param, n = 512, from = -3.5, 
     to = 0)
    } else {
      temp <- list(x = 3.5 * (seq(1, 512) - 1)/512, y = rep(0.001, 512))
    }
    obs.phi.densities.neg[, i] <- temp$y
    norm.factor <- sum(obs.phi.densities.neg[, i])
    obs.phi.densities.neg[, i] <- obs.phi.densities.neg[, i]/norm.factor
    
   }
   phi.density.mean.pos <- apply(phi.densities.pos, 1, mean)
   phi.density.mean.neg <- apply(phi.densities.neg, 1, mean)
   
   obs.phi.density.mean.pos <- apply(obs.phi.densities.pos, 1, mean)
   obs.phi.density.mean.neg <- apply(obs.phi.densities.neg, 1, mean)
   
   phi.density.median.pos <- apply(phi.densities.pos, 1, median)
   phi.density.median.neg <- apply(phi.densities.neg, 1, median)
   
   obs.phi.density.median.pos <- apply(obs.phi.densities.pos, 1, median)
   obs.phi.density.median.neg <- apply(obs.phi.densities.neg, 1, median)
   
   x <- c(x.coor.neg, x.coor.pos)
   x.plot.range <- range(x)
   y1 <- c(phi.density.mean.neg, phi.density.mean.pos)
   y2 <- c(obs.phi.density.mean.neg, obs.phi.density.mean.pos)
   y.plot.range <- c(-0.3 * max(c(y1, y2)), max(c(y1, y2)))
   
   print(c(y.plot.range, max(c(y1, y2)), max(y1), max(y2)))
   
   plot(x, y1, xlim = x.plot.range, ylim = 1.5 * y.plot.range, type = "l", 
    lwd = 2, col = 2, xlab = "NES", ylab = "P(NES)", main = "Global Observed and Null Densities (Area Normalized)")
   
   y1.point <- y1[seq(1, length(x), 2)]
   y2.point <- y2[seq(2, length(x), 2)]
   x1.point <- x[seq(1, length(x), 2)]
   x2.point <- x[seq(2, length(x), 2)]
   
   # for (i in 1:length(x1.point)) { lines(c(x1.point[i], x1.point[i]), c(0,
   # y1.point[i]), lwd = 3, cex = 0.9, col = colors()[555]) # shading } for (i in
   # 1:length(x2.point)) { lines(c(x2.point[i], x2.point[i]), c(0, y2.point[i]), lwd
   # = 3, cex = 0.9, col = colors()[29]) # shading }
   
   points(x, y1, type = "l", lwd = 2, col = colors()[555])
   points(x, y2, type = "l", lwd = 2, col = colors()[29])
   
   for (i in 1:Ng) {
    col <- ifelse(Obs.ES.norm[i] > 0, 2, 3)
    lines(c(Obs.ES.norm[i], Obs.ES.norm[i]), c(-0.2 * max(c(y1, y2)), 
      0), lwd = 1, lty = 1, col = 1)
   }
   leg.txt <- paste("Neg. ES: \"", phen2, " \" ", sep = "", collapse = "")
   text(x = x.plot.range[1], y = -0.25 * max(c(y1, y2)), adj = c(0, 1), 
    labels = leg.txt, cex = 0.9)
   leg.txt <- paste(" Pos. ES: \"", phen1, "\" ", sep = "", collapse = "")
   text(x = x.plot.range[2], y = -0.25 * max(c(y1, y2)), adj = c(1, 1), 
    labels = leg.txt, cex = 0.9)
   
   leg.txt <- c("Null Density", "Observed Density", "Observed NES values")
   c.vec <- c(colors()[555], colors()[29], 1)
   lty.vec <- c(1, 1, 1)
   lwd.vec <- c(2, 2, 2)
   legend(x = 0, y = 1.5 * y.plot.range[2], bty = "n", bg = "white", legend = leg.txt, 
    lty = lty.vec, lwd = lwd.vec, col = c.vec, cex = 0.9)
   if (gsea.type == "GSEA") {
    B <- A[obs.index, ]
    if (N > 300) {
      C <- rbind(B[1:100, ], rep(0, Ns), rep(0, Ns), B[(floor(N/2) - 
     50 + 1):(floor(N/2) + 50), ], rep(0, Ns), rep(0, Ns), B[(N - 
     100 + 1):N, ])
    }
    rm(B)
    GSEA.HeatMapPlot(V = C, col.labels = class.labels, col.classes = class.phen, 
      main = "Heat Map for Genes in Dataset")
   }
   # p-vals plot
   nom.p.vals <- p.vals[Obs.ES.index, 1]
   FWER.p.vals <- p.vals[Obs.ES.index, 2]
   plot.range <- 1.25 * range(NES)
   plot(NES, FDR.mean, ylim = c(0, 1), xlim = plot.range, col = 1, bg = 1, 
    type = "p", pch = 22, cex = 0.75, xlab = "NES", main = "p-values vs. NES", 
    ylab = "p-val/q-val")
   points(NES, nom.p.vals, type = "p", col = 2, bg = 2, pch = 22, cex = 0.75)
   points(NES, FWER.p.vals, type = "p", col = colors()[577], bg = colors()[577], 
    pch = 22, cex = 0.75)
   leg.txt <- c("Nominal p-value", "FWER p-value", "FDR q-value")
   c.vec <- c(2, colors()[577], 1)
   pch.vec <- c(22, 22, 22)
   legend(x = -0.5, y = 0.5, bty = "n", bg = "white", legend = leg.txt, 
    pch = pch.vec, col = c.vec, pt.bg = c.vec, cex = 0.9)
   lines(c(min(NES), max(NES)), c(nom.p.val.threshold, nom.p.val.threshold), 
    lwd = 1, lty = 2, col = 2)
   lines(c(min(NES), max(NES)), c(fwer.p.val.threshold, fwer.p.val.threshold), 
    lwd = 1, lty = 2, col = colors()[577])
   lines(c(min(NES), max(NES)), c(fdr.q.val.threshold, fdr.q.val.threshold), 
    lwd = 1, lty = 2, col = 1)
   
   dev.off()
   
  }  # if Ng > 1
 
 #----------------------------------------------------------------------------
 # Produce report for each gene set passing the nominal, FWER or FDR test or the
 # top topgs in each side
 
 if (topgs > floor(Ng/2)) {
  topgs <- floor(Ng/2)
 }
 
 for (i in 1:Ng) {
  if ((p.vals[i, 1] <= nom.p.val.threshold) || (p.vals[i, 2] <= fwer.p.val.threshold) || 
   (FDR.mean.sorted[i] <= fdr.q.val.threshold) || (is.element(i, c(Obs.ES.index[1:topgs], 
   Obs.ES.index[(Ng - topgs + 1):Ng])))) 
   {
    
    # produce report per gene set
    
    kk <- 1
    gene.number <- vector(length = size.G[i], mode = "character")
    gene.names <- vector(length = size.G[i], mode = "character")
    gene.symbols <- vector(length = size.G[i], mode = "character")
    gene.descs <- vector(length = size.G[i], mode = "character")
    gene.list.loc <- vector(length = size.G[i], mode = "numeric")
    core.enrichment <- vector(length = size.G[i], mode = "character")
    gene.rnk <- vector(length = size.G[i], mode = "numeric")
    gene.RES <- vector(length = size.G[i], mode = "numeric")
    rank.list <- seq(1, N)
    
    if (Obs.ES[i] >= 0) {
      set.k <- seq(1, N, 1)
      phen.tag <- phen1
      loc <- match(i, Obs.ES.index)
    } else {
      set.k <- seq(N, 1, -1)
      phen.tag <- phen2
      loc <- Ng - match(i, Obs.ES.index) + 1
    }
    
    for (k in set.k) {
      if (Obs.indicator[i, k] == 1) {
     gene.number[kk] <- kk
     # gene.names[kk] <- obs.gene.labels[k]
     gene.symbols[kk] <- obs.gene.symbols[k]
     gene.descs[kk] <- obs.gene.descs[k]
     gene.list.loc[kk] <- k
     gene.rnk[kk] <- signif(obs.rnk[k], digits = 3)
     gene.RES[kk] <- signif(Obs.RES[i, k], digits = 3)
     if (Obs.ES[i] >= 0) {
       core.enrichment[kk] <- ifelse(gene.list.loc[kk] <= Obs.arg.ES[i], 
      "YES", "NO")
     } else {
       core.enrichment[kk] <- ifelse(gene.list.loc[kk] > Obs.arg.ES[i], 
      "YES", "NO")
     }
     kk <- kk + 1
      }
    }
    
    gene.report <- data.frame(cbind(gene.number, gene.symbols, gene.descs, 
      gene.list.loc, gene.rnk, gene.RES, core.enrichment))
    if (gsea.type == "GSEA") {
      if (rank.metric == "S2N") {
     names(gene.report) <- c("#", "GENE SYMBOL", "DESC", "LIST LOC", 
       "S2N", "RES", "CORE_ENRICHMENT")
      }
      if (rank.metric == "ttest") {
     names(gene.report) <- c("#", "GENE SYMBOL", "DESC", "LIST LOC", 
       "TTest", "RES", "CORE_ENRICHMENT")
      }
      if (rank.metric == "seq") {
     names(gene.report) <- c("#", "GENE SYMBOL", "DESC", "LIST LOC", 
       "DESeq2 Log2(FC)", "RES", "CORE_ENRICHMENT")
      }
    } else if (gsea.type == "preranked") {
      names(gene.report) <- c("#", "GENE SYMBOL", "DESC", "LIST LOC", 
     "RNK", "RES", "CORE_ENRICHMENT")
    }
    
    # print(gene.report)
    
    if (output.directory != "") {
      
      filename <- paste(output.directory, doc.string, ".", gs.names[i], 
     ".report.", phen.tag, ".", loc, ".txt", sep = "", collapse = "")
      write.table(gene.report, file = filename, quote = FALSE, row.names = FALSE, 
     sep = "\t", na = "")
      
      gs.filename <- paste(output.directory, doc.string, ".", gs.names[i], 
     ".plot.", phen.tag, ".", loc, ".pdf", sep = "", collapse = "")
      pdf(file = gs.filename, height = 6, width = 14)
      
    }
    
    nf <- layout(matrix(c(1, 2, 3), 1, 3, byrow = T), 1, c(1, 1, 1), 
      TRUE)
    ind <- 1:N
    min.RES <- min(Obs.RES[i, ])
    max.RES <- max(Obs.RES[i, ])
    if (max.RES < 0.3) 
      max.RES <- 0.3
    if (min.RES > -0.3) 
      min.RES <- -0.3
    delta <- (max.RES - min.RES) * 0.5
    min.plot <- min.RES - 2 * delta
    max.plot <- max.RES
    max.corr <- max(obs.rnk)
    min.corr <- min(obs.rnk)
    Obs.correl.vector.norm <- (obs.rnk - min.corr)/(max.corr - min.corr) * 
      1.25 * delta + min.plot
    zero.corr.line <- (-min.corr/(max.corr - min.corr)) * 1.25 * delta + 
      min.plot
    col <- ifelse(Obs.ES[i] > 0, 2, 4)
    
    # Running enrichment plot
    
    sub.string <- paste("Number of genes: ", N, " (in list), ", size.G[i], 
      " (in gene set)", sep = "", collapse = "")
    
    main.string <- paste("Enrichment plot:", gs.names[i])
    plot(ind, Obs.RES[i, ], main = main.string, sub = sub.string, xlab = "Gene List Index", 
      ylab = "Running Enrichment Score (RES)", xlim = c(1, N), ylim = c(min.plot, 
     max.plot), type = "l", lwd = 2, cex = 1, col = col)
    for (j in seq(1, N, 20)) {
      lines(c(j, j), c(zero.corr.line, Obs.correl.vector.norm[j]), lwd = 1, 
     cex = 1, col = colors()[12])  # shading of correlation plot
    }
    lines(c(1, N), c(0, 0), lwd = 1, lty = 2, cex = 1, col = 1)  # zero RES line
    lines(c(Obs.arg.ES[i], Obs.arg.ES[i]), c(min.plot, max.plot), lwd = 1, 
      lty = 3, cex = 1, col = col)  # max enrichment vertical line
    for (j in 1:N) {
      if (Obs.indicator[i, j] == 1) {
     lines(c(j, j), c(min.plot + 1.25 * delta, min.plot + 1.75 * delta), 
       lwd = 1, lty = 1, cex = 1, col = 1)  # enrichment tags
      }
    }
    lines(ind, Obs.correl.vector.norm, type = "l", lwd = 1, cex = 1, 
      col = 1)
    lines(c(1, N), c(zero.corr.line, zero.corr.line), lwd = 1, lty = 1, 
      cex = 1, col = 1)  # zero correlation horizontal line
    temp <- order(abs(obs.rnk), decreasing = T)
    arg.correl <- temp[N]
    lines(c(arg.correl, arg.correl), c(min.plot, max.plot), lwd = 1, 
      lty = 3, cex = 1, col = 3)  # zero crossing correlation vertical line
    
    leg.txt <- paste("\"", phen1, "\" ", sep = "", collapse = "")
    text(x = 1, y = min.plot, adj = c(0, 0), labels = leg.txt, cex = 1)
    
    leg.txt <- paste("\"", phen2, "\" ", sep = "", collapse = "")
    text(x = N, y = min.plot, adj = c(1, 0), labels = leg.txt, cex = 1)
    
    adjx <- ifelse(Obs.ES[i] > 0, 0, 1)
    
    leg.txt <- paste("Peak at ", Obs.arg.ES[i], sep = "", collapse = "")
    text(x = Obs.arg.ES[i], y = min.plot + 1.8 * delta, adj = c(adjx, 
      0), labels = leg.txt, cex = 1)
    
    leg.txt <- paste("Zero crossing at ", arg.correl, sep = "", collapse = "")
    text(x = arg.correl, y = min.plot + 1.95 * delta, adj = c(adjx, 0), 
      labels = leg.txt, cex = 1)
    
    # nominal p-val histogram
    
    sub.string <- paste("ES =", signif(Obs.ES[i], digits = 3), " NES =", 
      signif(Obs.ES.norm[i], digits = 3), "Nom. p-val=", signif(p.vals[i, 
     1], digits = 3), "FWER=", signif(p.vals[i, 2], digits = 3), "FDR=", 
      signif(FDR.mean.sorted[i], digits = 3))
    temp <- density(phi[i, ], adjust = adjust.param)
    x.plot.range <- range(temp$x)
    y.plot.range <- c(-0.125 * max(temp$y), 1.5 * max(temp$y))
    plot(temp$x, temp$y, type = "l", sub = sub.string, xlim = x.plot.range, 
      ylim = y.plot.range, lwd = 2, col = 2, main = "Gene Set Null Distribution", 
      xlab = "ES", ylab = "P(ES)")
    x.loc <- which.min(abs(temp$x - Obs.ES[i]))
    lines(c(Obs.ES[i], Obs.ES[i]), c(0, temp$y[x.loc]), lwd = 2, lty = 1, 
      cex = 1, col = 1)
    lines(x.plot.range, c(0, 0), lwd = 1, lty = 1, cex = 1, col = 1)
    
    leg.txt <- c("Gene Set Null Density", "Observed Gene Set ES value")
    c.vec <- c(2, 1)
    lty.vec <- c(1, 1)
    lwd.vec <- c(2, 2)
    legend(x = -0.2, y = y.plot.range[2], bty = "n", bg = "white", legend = leg.txt, 
      lty = lty.vec, lwd = lwd.vec, col = c.vec, cex = 1)
    
    leg.txt <- paste("Neg. ES \"", phen2, "\" ", sep = "", collapse = "")
    text(x = x.plot.range[1], y = -0.1 * max(temp$y), adj = c(0, 0), 
      labels = leg.txt, cex = 1)
    leg.txt <- paste(" Pos. ES: \"", phen1, "\" ", sep = "", collapse = "")
    text(x = x.plot.range[2], y = -0.1 * max(temp$y), adj = c(1, 0), 
      labels = leg.txt, cex = 1)
    
    # create pinkogram for each gene set
    
    kk <- 1
    if (gsea.type == "GSEA") {
      pinko <- matrix(0, nrow = size.G[i], ncol = cols)
      pinko.gene.names <- vector(length = size.G[i], mode = "character")
      for (k in 1:rows) {
     if (Obs.indicator[i, k] == 1) {
       pinko[kk, ] <- A[obs.index[k], ]
       pinko.gene.names[kk] <- obs.gene.symbols[k]
       kk <- kk + 1
     }
      }
      GSEA.HeatMapPlot(V = pinko, row.names = pinko.gene.names, col.labels = class.labels, 
     col.classes = class.phen, col.names = sample.names, main = " Heat Map for Genes in Gene Set", 
     xlab = " ", ylab = " ")
    }
    dev.off()
    
   }  # if p.vals thres
  
 }  # loop over gene sets
 
 return(list(report1 = report.phen1, report2 = report.phen2))
 
}
