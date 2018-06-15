#' utility fn for getComplexity and preseqR 
#'
#' @param x   the thing to count features in (usually a GRanges)
#'
#' @return    an array of feature frequencies
#' 
#' @import    GenomicRanges
#' 
#' @export
getFeatureCounts <- function(x) {

  if (!is(x, 'GRanges')) x <- as(x, 'GRanges')
  ends <- paste(start(x), end(x), sep=':')
  freqs <- data.frame(seen=as(table(ends), 'matrix')[,1])
  freqs[ rev(order(freqs$seen)), , drop=FALSE] ## descending order 

}
