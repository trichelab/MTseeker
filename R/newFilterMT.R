#' Filter SummarizedExperiment, DataFrame, MVRanges, or MVRangesList using user specified constraints
#'
#' Griffin et al. (Genetics in Medicine 2014) recommends 20x coverage for mtDNA
#' sequencing to have comparable error rates to Sanger sequencing.  By default, 
#' that is the cutoff applied here to ensure halfway decent variant annotation.
#' 
#' Triska (Cancer Res, in revision) suggests a small number of masked regions 
#' where homopolymers can be a problem; these are avoided if fpFilter
#' 
#' The NuMT filtration step (Ju, in eLife 2014, suggests a variant allele cutoff
#' of 0.03 to avoid false positive calls from nuclear-mitochondrial translocated
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
#'
#' @param vars              variants, can be MVRanges[List] or DataFrame/SummarizedExperiment with colData()$`mtCovg`
#' @param minTotalDepth     minimum number of total reads (ref + alt reads) (20)
#' @param minAltDepth       minimum number of reads supporting a mutation needed (2)
#' @param minVAF            minimum VAF to be considered (0.90)
#' @param minCovg           minimum total depth (20, cf. Griffin, Genetics in Medicine 2014)
#' @param fpFilter          apply Triska's homopolymer false positive filter? Only applicable for rCRS (FALSE)
#' @param NuMT              apply the 0.03 VAF NuMT filter from Ju (GR 2015)? (FALSE)
#' @param verbose           prints the filters applied (FALSE)
#'
#' @return            a filtered SE, data.frame, MVRanges, or MVRangesList 
#' 
#' @import SummarizedExperiment
#' @import S4Vectors
#'
#' @examples
#' filterMT(data.frame(sample="foo", mtCovg=1000))
#'
#' @export
newFilterMT <- function(vars, minTotalDepth=20, minAltDepth=2, minVAF=0.90,
                         minCovg=20, fpFilter=TRUE, NuMT=TRUE, verbose=FALSE) {
  
  if (is(vars, "MVRangesList")) {
    
    # Since this is a very customizable function
    # I thought it might be a good idea to print what the settings have been set to
    if (verbose) {
      
      if (NuMT) message("Discarding calls with VAF < 3% to avoid NuMT contamination...")
      if (fpFilter) message("Filtering out common false positive regions in rCRS...") 
      message("Filtering out samples with < ", minCovg, "x mean read coverage...")
      message("Discarding calls with a total depth < ", minTotalDepth, " and alt depth < ", minAltDepth, "...")
    }

    mvr <- lapply(vars, newFilterMT, minTotalDepth=minTotalDepth, minAltDepth=minAltDepth, minVAF=minVAF,
                                      minCovg=minCovg, fpFilter=fpFilter, NuMT=NuMT, verbose=verbose)
    return(MVRangesList(mvr))
  } 
  
  if (is(vars, "MVRanges")) { 

    # If a sample has on average, less than minCovg coverage (default is 20)
    # Then it is a low quality sample that should be discarded
    if (genomeCoverage(vars) < minCovg) {
      return(vars[0])
    }
    
    # More complicated version for human genome
    if (genome(vars) == "rCRS") {
      
      return(.rCRSfilter(vars, fp=fpFilter, NuMT=NuMT, 
                         minTotalDepth=minTotalDepth, minAltDepth=minAltDepth, minVAF=minVAF))
    }
    
    else {
      
      # Only keep variants with a minimum altDept and totalDepth as specified by the user
      vars <- subset(vars, !is.na(altDepth(vars)) & altDepth(vars) >= minAltDepth)
      vars <- subset(vars, !is.na(totalDepth(vars)) & totalDepth(vars) >= minTotalDepth)
      vars <- subset(vars, VAF >= minVAF)
      
      # Nuclear contamination
      if (NuMT) vars <- subset(vars, VAF >= 0.03)
      
      return(vars)
    }
  } #MVRanges
  
  # Accepts data frames as well
  else if (is(vars, "DataFrame") | is(vars, "data.frame")) {
    stop("Gotta fix for dataframes")
  }
  
  else {
    stop("filterMT operates only on MVRangesLists, MVRanges, and data frames.")
  }
  
}

# helper fn
# If the input is a data frame
.variantsDF <- function(DFSE, fpFilter=FALSE, NuMT=TRUE, verbose=FALSE) {
  
  if (is(DFSE, "DataFrame") | is(DFSE, "data.frame")) {
    stopifnot("mtCovg" %in% names(DFSE))
  } else if (is(DFSE, "SummarizedExperiment") & 
             !"mtCovg" %in% names(colData(DFSE))) {
    stop("filterMT requires colData named `mtCovg`.")
  } else if (!"mtCovg" %in% names(DFSE)) {
    stop("filterMT requires a column named `mtCovg`.")
  }
  
  if(verbose) message("Filtering out samples with < ", minCovg, "x mean read coverage...")
  if (is(DFSE, "SummarizedExperiment")) {
    res <- DFSE[, which(DFSE$mtCovg >= minCovg)]
  } else { 
    res <- subset(DFSE, mtCovg >= minCovg)
  }
  
  if (fpFilter) {
    if (verbose) message("Filtering out common false positive regions...") 
    data(fpFilter_Triska, package="MTseeker") 
    res <- subsetByOverlaps(res, subset(gaps(fpFilter_Triska), strand=="*"))
  }
  
  if (NuMT) {
    if (is(res, "SummarizedExperiment")) {
      if (verbose) message("Discarding calls with VAF < 3% to avoid NuMT contamination...")
      res <- subset(res, res$VAF >= 0.03)
    }  
  }
  
  return(res)
}

# helper fn
# More complicated filter for rCRS ref genome
.rCRSfilter <- function(vars, fp=TRUE, NuMT=TRUE, minTotalDepth=20, minAltDepth=2, minVAF=0.90) {

  # False positive filter
  if (fp) {
    data(fpFilter_Triska, package="MTseeker")
    fpFilter <- subset(gaps(fpFilter_Triska), strand == "*")
  }
  
  # a nonfilter -- keep anything and everything on chrM
  else {
    fpFilter <- GRanges("chrM", IRanges(start=1, end=16569), strand="*")
  }
  
  # Only keep variants with an altDept of at least minAltDepth
  vars <- subset(vars, !is.na(altDepth(vars)) & altDepth(vars) >= minAltDepth)
  vars <- subset(vars, !is.na(totalDepth(vars)) & totalDepth(vars) >= minTotalDepth)
  vars <- subset(vars, VAF >= minVAF)
  
  # Nuclear contamination
  if (NuMT) {
    vars <- subset(subsetByOverlaps(vars, fpFilter), VAF >= 0.03)
  } else {
    vars <- subsetByOverlaps(vars, fpFilter)
  }
  
  #if ("PASS" %in% names(mcols(vars))) vars <- subset(vars, PASS)

  # If it is empty, return the empty MVRanges
  if (length(vars) == 0) return(vars)
  
  names(vars) <- MTHGVS(vars)

  return(vars)
}
