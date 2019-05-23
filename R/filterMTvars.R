#' sanitize PASSing mitochondrial variant calls to a moderate degree
#'
#' @param vars 	an MVRanges or MVRangesList (will be unlisted and relisted)
#' @param fp 	  use false positive filter[s]? (TRUE: use Triska fpFilter)
#' @param NuMT	variants with VAF < [this number] will be presumed NuMTs (0.03)
#' @param covg	minimum median read coverage across chrM to be considered (20) 
#' @param depth minimum altDepth for a variant to be retained (2) 
#'
#' @return	a filtered set of variants 
#' 
#' @import 	GenomicRanges
#'
#' @examples
#' library(MTseekerData) 
#' filterMTvars(RONKSvariants$RO_1)
#'
#' @export 
filterMTvars <- function(vars, fp=TRUE, NuMT=0.03, covg=20, depth=2) {

  if (fp) { 
    data(fpFilter_Triska, package="MTseeker")
    fpFilter <- subset(gaps(fpFilter_Triska), strand == "*")
  } else { 
    # a nonfilter -- keep anything and everything on chrM
    fpFilter <- GRanges("chrM", IRanges(start=1, end=16569), strand="*")
  } 

  if (is(vars, "MVRanges")) {
    if (!"VAF" %in% names(mcols(vars))) {
      # add VAF column for NuMT filtering
      vars$VAF <- altDepth(vars) / totalDepth(vars)
    }
    if (!"PASS" %in% names(mcols(vars))) {
      # mark all vars as PASS if unflagged (!)
      # currently, this happens with pileupMT()
      vars$PASS <- TRUE # this is sketchy AF
    }
    subset(subsetByOverlaps(vars, fpFilter), 
           VAF >= NuMT & PASS & altDepth(vars) >= depth)
  } else if (is(vars, "MVRangesList")) {
    MVRangesList(lapply(vars[genomeCoverage(vars)>=covg], 
                        filterMTvars, fp=fp, NuMT=NuMT))
  } else { 
    stop("This function is only meant for MVRanges and MVRangesList objects.")
  }

}
