#' Inject (one or more) variants against rCRS.
#' 
#' FIXME: this function could most likely be orders of magnitude faster.
#' FIXME: this ONLY considers variants injected against rCRS, not RSRS or hg19. 
#' 
#' @param mvr       An MVRanges, usually from callMT, often subsetted
#' @param gr        A GRanges, usually of protein-coding regions (the default)
#' @param xlate     Attempt to translate codon(s) affected by variant(s)? (TRUE)
#' @param canon     Minimum VAF to treat variants as canonical by subject (0.99)
#' @param refX      Reference depth below which variant is deemed canonical (1)
#' @param altX      Alternative depth above which variants deemed canonical (1)
#' 
#' @return          The GRanges, with ref/var DNA and AA and 
#' 
#' @import GenomicRanges 
#'
#' @export
injectMtVariants <- function(mvr, gr=NULL, xlate=TRUE, 
                             canon=.99, refX=1, altX=1) {

  # rCRS only, for the time being 
  stopifnot(unique(genome(mvr)) == "rCRS")

  # get mtGenes if needed 
  if (is.null(gr)) gr <- genes(mvr)  
  stopifnot(unique(genome(gr)) == "rCRS")

  # subset the variants to those that overlap the target GRanges and are canon
  mvr <- subset(locateVariants(subsetByOverlaps(mvr, gr, type="within")),
                VAF >= canon & refDepth < refX & altDepth > altX )

  # mitochondrial genomic sequence
  data(rCRSeq, package="MTseeker")
  gr$refDNA <- getSeq(rCRSeq, gr)
  altSeq <- DNAStringSet(replaceAt(rCRSeq[[1]], ranges(mvr), alt(mvr)))
  names(altSeq) <- names(rCRSeq)
  gr$varDNA <- getSeq(altSeq, gr)

  if (xlate) {

    # use MT_CODE to translate results
    MT_CODE <- getGeneticCode("SGC1")
    gr$refSeq <- getSeq(rCRSeq, gr)
    gr$refAA <- suppressWarnings(translate(gr$refSeq, MT_CODE))

    # this is tricky!
    for (g in seq_along(gr)) {
      submvr <- subsetByOverlaps(mvr, gr[g])
      gr[g]$varSeq <- replaceAt(gr[g]$refSeq, 
                                IRanges(submvr$localStart, submvr$localEnd),
                                alt(submvr))
    } 

    gr$varAA <- suppressWarnings(translate(gr$varDNA, MT_CODE))


    # which codon(s) have been perturbed?
    alignments <- pairwiseAlignment(gr$refAA, 
                                    gr$varAA,
                                    substitutionMatrix = "BLOSUM50",
                                    gapOpening = 0, gapExtension = 8)


      
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
