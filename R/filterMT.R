#' Filter SummarizedExperiment or DataFrame values based on its $mtCovg column
#'
#' Griffin et al. (Genetics in Medicine 2014) recommends 20x coverage for mtDNA
#' sequencing to have comparable error rates to Sanger sequencing.  By default, 
#' that is the cutoff applied here to ensure halfway decent variant annotation.
#' 
#' Triska (Cancer Res, in revision) suggests a small number of masked regions 
#' where homopolymers can be a problem; these are avoided if fpFilter
#' 
#' The NuMT filtration step (Ju, in Genome Research 2015, suggests a 2.5% VAF 
#' cutoff to avoid false positive calls from nuclear-mitochondrial translocated
#' or 'NuMT' fragments) is also a useful tool to cut down on nonsensical calls,
#' although it may be important to use caution as low heteroplasmy can also 
#' resolve into apparent near-homoplasmy at the single-cell level, at least in
#' our (TJT & co) experience.
#' 
#' As a consequence of the Wild West nature for published methods of high-
#' throughput mitochondrial sequence variant analysis at the time of writing 
#' (2018), the default for this function is to filter on coverage only, as 
#' the user is expected to determine what additional filters to apply. We
#' could envision changing these defaults down the road as standards congeal.
#'
#' @param DFSE        a DataFrame/SummarizedExperiment with colData()$`mtCovg`
#' @param minCovg     minimum covg (20, cf. Griffin, Genetics in Medicine 2014)
#' @param fpFilter    apply Triska's homopolymer false positive filter? (FALSE)
#' @param NuMT        apply the 2.5% VAF NuMT filter from Ju (GR 2015)? (FALSE)
#'
#' @import SummarizedExperiment
#' @import S4Vectors
#'
#' @export
filterMT <- function(DFSE, minCovg=20, fpFilter=FALSE, NuMT=FALSE) {

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
    res <- DFSE[, which(DFSE$mtCovg >= minCovg)]
  } else { 
    res <- subset(DFSE, mtCovg >= minCovg)
  }

  if (fpFilter) {
    message("Filtering out common false positive regions...") 
    data(fpFilter_Triska, package="MTseeker") 
    res <- subsetByOverlaps(res, subset(gaps(fpFilter_Triska), strand=="*"))
  }

  if (NuMT) {
    if (is(res, "SummarizedExperiment")) {
      message("Discarding calls with VAF < 2.5% to avoid NuMT contamination...")
      res <- subset(res, res$VAF >= 0.025)
    }  
  }

  return(res)

}
