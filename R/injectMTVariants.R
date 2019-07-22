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
      
      if (length(submvr)) {
        
        ### FIX: Not sure what the actual codon end should be for the second gene
        ### FIX: For now this only handles deletions
        
        # Check to see if the variant is in an overlapping gene
        # Assumes the entire second gene will not be deleted
        if (!is.na(submvr$overlapGene)) {
          
          # Workaround 
          # Not fixed for insertions
          if (grepl("ins", names(submvr))) {
            
            orig <- extractAt(gr[g]$refAA, IRanges(submvr$startCodon,submvr$startCodon))
            altd <- extractAt(gr[g]$varAA, IRanges(submvr$startCodon,submvr$startCodon))
            gr[g]$consequences <- "overlap"
            warning("Overlapping genes with insertion in: ", print(names(submvr)))
            next
          }
          
          # Names of overlapping gene
          affectedGenes <- unlist(strsplit(mcols(submvr)$overlapGene, split=","))
          
          # For the first gene
          if (g == affectedGenes[1]) {

            # If the variant is located within the range of start and end of the gene
            # then it can continue as normal
            if ( (start(submvr) >= start(gr[g])) && (end(submvr) <= end(gr[g])) ) {
              # This is fine
            }
            
            # If the variant starts in one gene, but then goes into the next
            # then artificially makes a new ending for the variant within the first gene
            if (end(submvr) > end(gr[g])) {
              
              # Create an artifical end for the variant
              submvr$localEnd <- width(gr[g])
              submvr$endCodon <- width(gr[g]$refAA)
              
            }
          } # 1st gene
          
          # For the second gene
          # In this next iteration, all information about the variant is reset
          if (g == affectedGenes[2]) {
            
            # If the deletion is completely encompassed within this second gene
            if ( (start(submvr) >= start(gr[g])) && (end(submvr) <= end(gr[g])) ) {
              
              # Have to change the start and end locations to reflect local values
              submvr$localStart <- start(submvr) - start(gr[g])
              submvr$localEnd <- submvr$localStart + (width(submvr) - 1)
              
            }
            
            # If part of the deletion takes place in the first gene
            # then continues the deletion in the beginning of this second gene
            # Must alter the sequence that is replaced
            if (start(submvr) < start(gr[g])) {
              
              oldWidth <- width(submvr)
              submvr$localStart <- 1
              submvr$localEnd <- submvr$localStart + (oldWidth - 1) 
              
            }
            
            # Need to figure out the startCodon and endCodon
            # For now, if the sequence is evenly divided into 3
            # assume the first bp of the gene is the start bp
            if ( (width(gr[g]$refSeq) %% 3 == 0) && (width(gr[g]$refSeq) / 3 == width(gr[g]$refAA)) ){
              
              submvr$startCodon <- floor(submvr$localStart / 3)
              submvr$endCodon <- floor(submvr$localEnd / 3)
            }
            
            else {
              message("Need new method to determine start and end codons")
              warning("Setting codons potentially incorrectly for ", names(submvr))
              submvr$startCodon <- floor(submvr$localStart / 3)
              submvr$endCodon <- floor(submvr$localEnd / 3)
            }
            
            
          } # 2nd gene
          
        } # !overlappingGene
      }
      
      #related to the off by 1 error below
      #this should catch indels in codon 1
      if (length(submvr)) {
        if (submvr$localStart == 0 & submvr$localEnd > 0) {
          submvr$localStart <- 1
          submvr$localEnd <- submvr$localEnd + 1
        }
        if (submvr$localStart == 0 & submvr$localEnd == 0) {
          submvr$localStart <- 1
          submvr$localEnd <- 1
        }
      }
      subir <- IRanges(submvr$localStart, submvr$localEnd)
      gr[g]$varSeq <- replaceAt(gr[g]$refSeq, subir, alt(submvr))
      gr[g]$varAA <- suppressWarnings(translate(gr[g]$varSeq, MT_CODE))
      # FIXME: probably different in terms of end codon from orig (esp for AAfs)
      #bizarre checks for if a bp change occurs in the first codon at bp 1
      #this actually looks like an off by 1 error when the 1st codon is impacted
      if (length(submvr)) {
        if (submvr$startCodon == 0) {
          submvr$startCodon <- 1
        }
      }
      if (length(submvr)) {
        if (submvr$endCodon == 0) {
          submvr$endCodon <- 1
        }
      }
      
      # AA sequence from the reference
      orig <- extractAt(gr[g]$refAA, IRanges(submvr$startCodon,submvr$endCodon))
      
      # If there is a deletion at the end of a gene
      # Must go here otherwise the endCodon value changes and you get the incorrect reference AA
      if (length(submvr)) {
        if (submvr$endCodon > width(gr[g]$varAA)) {
          
          submvr$endCodon <- width(gr[g]$varAA)
        }
      }
      
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
