#' Inject (one or more) variants against rCRS.
#' 
#' FIXME: this function could most likely be orders of magnitude faster.
#' FIXME: this ONLY considers variants injected against rCRS, not RSRS or hg19. 
#' 
#' @param mvr       An MVRanges, usually from callMT, often subsetted
#' @param gr        A GRanges, usually of protein-coding regions (the default)
#' @param aa        Attempt to translate codon(s) affected by variant(s)? (TRUE)
#' @param canon     Minimum VAF to treat variants as canonical by subject (0.99)
#' @param refX      Reference depth below which variant is deemed canonical (1)
#' @param altX      Alternative depth above which variants deemed canonical (1)
#' 
#' @return          The GRanges, with ref/var DNA and AA and 
#' 
#' @import GenomicRanges 
#' @importFrom utils data
#'
#' @examples
#' library(MTseekerData)
#' RO_2 <- RONKSvariants[["RO_2"]]
#' injectMTVariants(RO_2)
#'
#' @export
injectMTVariants <- function(mvr, gr=NULL, aa=TRUE, canon=.99, refX=1, altX=1) {

  # rCRS only, for the time being 
  stopifnot(unique(genome(mvr)) == "rCRS")

  # get mtGenes if needed 
  if (is.null(gr)) gr <- genes(mvr)
  stopifnot(unique(genome(gr)) == "rCRS")

  # subset the variants to those that overlap the target GRanges and are canon
  mvr <- subset(locateVariants(subsetByOverlaps(mvr, gr, type="within")),
                VAF >= canon & refDepth < refX & altDepth > altX )

  # mitochondrial genomic sequence
  # FIXME: may not want to do this
  # FIXME: else may want to filter
  data(rCRSeq, package="MTseeker")
  gr$refDNA <- getSeq(rCRSeq, gr)
  altSeq <- DNAStringSet(replaceAt(rCRSeq[[1]], ranges(mvr), alt(mvr)))
  names(altSeq) <- names(rCRSeq)
  gr$varDNA <- getSeq(altSeq, gr)

  if (aa) {

    # use MT_CODE to translate results
    MT_CODE <- getGeneticCode("SGC1")
    gr$refSeq <- getSeq(rCRSeq, gr)
    gr$varSeq <- gr$refSeq 
    gr$refAA <- suppressWarnings(translate(gr$refSeq, MT_CODE))
    gr$varAA <- gr$refAA
    gr$consequences <- NA_character_

    # this is a bit tricky 
    for (g in names(gr)) {
      submvr <- subsetByOverlaps(mvr, gr[g])
      subir <- IRanges(submvr$localStart, submvr$localEnd)
      gr[g]$varSeq <- replaceAt(gr[g]$refSeq, subir, alt(submvr))
      gr[g]$varAA <- suppressWarnings(translate(gr[g]$varSeq, MT_CODE))
      # FIXME: probably different in terms of end codon from orig (esp for AAfs)
      #bizarre checks for if a bp change occurs in the first codon at bp 1
      # if (length(submvr)) {
      #   if (submvr$startCodon == 0) {
      #     submvr$startCodon <- 1
      #   }
      # }
      # if (length(submvr)) {
      #   if (submvr$endCodon == 0) {
      #     submvr$endCodon <- 1
      #   }
      # }
      orig <- extractAt(gr[g]$refAA, IRanges(submvr$startCodon,submvr$endCodon))
      altd <- extractAt(gr[g]$varAA, IRanges(submvr$startCodon,submvr$endCodon))
      gr[g]$consequences <- .flattenConsequences(orig, altd, submvr$startCodon)

    } 

    # for later
    if (FALSE) {
      alignments <- pairwiseAlignment(gr$refAA, gr$varAA,
                                      substitutionMatrix = "BLOSUM50",
                                      gapOpening = 0, gapExtension = 8)
    }
      
  }

  # done
  return(gr)
}

# helper function
.flattenConsequences <- function(orig, altd, startCodons) {

  .flat <- function(x) vapply(x[[1]], as.character, "N") 
  .prettify <- function(x) paste0(gsub(" ", "", x), collapse="")
  asdf <- DataFrame(refAA=.flat(orig), pos=startCodons, varAA=.flat(altd))
  csqs <- apply(subset(asdf, asdf$refAA != asdf$varAA), 1, .prettify)
  paste(csqs, collapse=",")

}
