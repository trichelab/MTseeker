#' Inject (one or more) variants against rCRS.
#' 
#' FIXME: this function could most likely be orders of magnitude faster.
#' FIXME: this ONLY considers variants injected against rCRS, not RSRS or hg19. 
#' 
#' @param gr        A GRanges, usually of protein-coding regions
#' @param mvr       An MVRanges, usually from callMT
#' @param xlate     Attempt to translate codon(s) affected by variant(s)? (TRUE)
#' 
#' @return          The same GRanges, but with ref/var DNA and AA from rCRS. 
#' 
#' @import GenomicRanges 
#'
#' @export
injectMtVariants <- function(gr, mvr, xlate=TRUE) {

  # rCRS only, for the time being 
  stopifnot(unique(genome(gr)) == "rCRS")
  stopifnot(unique(genome(mvr)) == "rCRS")
  
  # subset the variants to those that overlap the target GRanges
  mvr <- locateVariants(subsetByOverlaps(mvr, gr)) 

  # mitochondrial genomic sequence
  data(rCRSeq, package="MTseeker")
  gr$refDNA <- getSeq(rCRSeq, gr)
  altSeq <- DNAStringSet(replaceLetterAt(rCRSeq[[1]], start(mvr), alt(mvr)))
  names(altSeq) <- names(rCRSeq)
  gr$varDNA <- getSeq(altSeq, gr)

  if (xlate) {

    # which codon(s) have been perturbed?
    # split codons out to do this:
    stop("Not quite finished") 
      
    # use MT_CODE to translate results
    MT_CODE <- getGeneticCode("SGC1")
    gr$refAA <- suppressWarnings(translate(gr$refCodon, MT_CODE))
    gr$varAA <- suppressWarnings(translate(gr$varCodon, MT_CODE))

  }

  # done
  return(gr)
}

# helper function
.injectVars <- function(mvr) {
  data(rCRSeq, package="MTseeker")
  replaceAt(rCRSeq, ranges(mvr), alt(mvr))
}

# helper function 
.getCodons <- function(gr, mvr) { 
  data(rCRSeq, package="MTseeker")
  gr$firstCodon <- (start(mvr) - start(gr)) %/% 3
  gr$lastCodon <- (end(mvr) - start(gr)) %/% 3
  return(gr)
}
