#' Filter SummarizedExperiment or DataFrame values based on its $mtCovg column
#'
#' Griffin et al. (Genetics in Medicine 2014) recommends 20x coverage for mtDNA
#' sequencing to have comparable error rates to Sanger sequencing.  By default, 
#' that is the cutoff applied here to ensure halfway decent variant annotation.
#' 
#' @param DFSE        a DataFrame/SummarizedExperiment with colData() $mtCovg
#' @param minCovg     minimum covg (20, cf. Griffin, Genetics in Medicine 2014)
#'
#' @import S4Vectors
#' @import SummarizedExperiment
#'
#' @export
filterMT <- function(DFSE, minCovg=20) {

  if (!class(DFSE) %in% 
      c("data.frame","DataFrame",
        "SummarizedExperiment","RangedSummarizedExperiment")) {
    stop("filterMT operates on a [Ranged]SummarizedExperiment or DataFrame.")
  } else {
    if (is(DFSE, "DataFrame") | is(DFSE, "data.frame")) {
      stopifnot("mtCovg" %in% names(DFSE))
    } else if (is(DFSE, "SummarizedExperiment") & 
               !"mtCovg" %in% names(colData(DFSE))) {
      stop("filterMT requires colData named `mtCovg`.")
    } else if (!"mtCovg" %in% names(DFSE)) {
      stop("filterMT requires a column named `mtCovg`.")
    }
  }

  message("Filtering out samples with < ", minCovg, "x mean read coverage...")
  if (is(DFSE, "SummarizedExperiment")) {
    return(DFSE[, which(DFSE$mtCovg >= minCovg)])
  } else { 
    return(subset(DFSE, mtCovg >= minCovg))
  }

}
