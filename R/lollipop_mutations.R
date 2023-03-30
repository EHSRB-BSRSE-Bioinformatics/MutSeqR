#' Plot mutations in lollipop plot
#' 
#' Uses the trackViewer package to plot mutations in a lollipop plot in specific
#' regions as defined by the user input.
#' 
#' @param species One of "human" or "mouse"
#' @param mutations A GRanges object with mutation data
#' @param ... Additional arguments to trackViewer::lolliplot (e.g.,
#' `ranges = GRanges("chr1", IRanges(104, 109))` )
#' @importFrom trackViewer lolliplot
#' @importFrom GenomicRanges GRanges
#' @export
lollipop_mutations <- function(species = "human",
                               mutations = mutation_data,
                               ...
                               ) {
  if (!class(mutations) == "GRanges") { stop("Please supply a GRanges 
                                                object as the `mutations`
                                                parameter.")}
  features.gr <- get_region_seqs(species)
  SNP.gr <- mutations[elementMetadata(mutations)[["variation_type"]] == "snv"]
  SNP.gr$names <- SNP.gr$normalized_context_with_mutation
  SNP.gr$color <- c("black")
  #SNP.gr$border <- c("black")
  #SNP.gr$score <- SNP.gr$alt_depth
  #SNP.gr$score <- sample.int(5, length(SNP.gr), replace = TRUE)
  #SNP.gr.1 <- plyranges::filter_by_overlaps(SNP.gr, features.gr.1)
  features.gr.1 <- features.gr[1]
  lol <- trackViewer::lolliplot(SNP.gr, features.gr.1) #, ...)
  return(lol)
}

# features <- GRanges("chr1", IRanges(c(1, 501, 1001), 
#                                     width=c(120, 400, 405),
#                                     names=paste0("block", 1:3)),
#                     fill = c("#FF8833", "#51C6E6", "#DFA32D"),
#                     height = c(0.02, 0.05, 0.08))
# SNP <- c(10, 100, 105, 108, 400, 410, 420, 600, 700, 805, 840, 1400, 1402)
# sample.gr <- GRanges("chr1", IRanges(SNP, width=1, names=paste0("snp", SNP)),
#                      color = sample.int(6, length(SNP), replace=TRUE),
#                      score = sample.int(5, length(SNP), replace = TRUE))
# trackViewer::lolliplot(sample.gr, features)
