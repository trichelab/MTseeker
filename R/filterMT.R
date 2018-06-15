#' Filter a RangedSummarizedExperiment's columns based on its $mtCovg colData 
#'
#' Griffin et al. (Genetics in Medicine 2014) recommends 20x coverage for mtDNA
#' sequencing to have comparable error rates to Sanger sequencing.  By default, 
#' that is the cutoff applied here to ensure halfway decent variant annotation.
#' 
#' @param SE          a RangedSummarizedExperiment with colData() named mtCovg
#' @param minCovg     minimum covg (20, cf. Griffin, Genetics in Medicine 2014)
#'
#' @import SummarizedExperiment
#'
#' @export
filterMT <- function(SE, minCovg=20) {

  if (!is(SE, "RangedSummarizedExperiment")) {
    stop("filterMT operates on a RangedSummarizedExperiment.")
  } else if (! "mtCovg" %in% names(colData(SE)) ) { 
    stop("filterMT requires a RangedSummarizedExperiment with colData $mtCovg.")
  }

  message("Filtering out samples with < ", minCovg, "x mean read coverage...")
  return(SE[, which(SE$mtCovg >= minCovg)])

}
