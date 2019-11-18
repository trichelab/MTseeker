#' Inject (one or more) variants against rCRS.
#' 
#' FIXME: this function could most likely be orders of magnitude faster.
#' FIXME: this ONLY considers variants injected against rCRS, not RSRS or hg19. 
#' FIXME: this function really needs mouse support as well.
#' 
#' @param mvr       An MVRanges, usually from pileupMT, often subsetted
#' @param gr        A GRanges, usually of protein-coding regions (the default)
#' @param coding    TRUE to only look at coding regions
#' 
#' @return          The GRanges, with ref/var DNA and AA and 
#' 
#' @import GenomicRanges 
#' @importFrom utils data
#'
#' @examples
#' library(MTseekerData)
#' injectMTVariants(RONKSvariants[["RO_2"]]) 
#' 
#' @export
injectMTVariants <- function(mvr, coding=TRUE, gr=NULL) {

  # rCRS only, for the time being 
  # FIXME: add mouse, at the very least
  stopifnot(unique(genome(mvr)) == "rCRS")
  
  # get mtGenes if needed 
  if (is.null(gr)) gr <- .getAnno(genome(mvr))
  stopifnot(unique(genome(gr)) == "rCRS")
  
  # Assuming you have filtered by now
  #mvr <- newFilterMT(mvr, minTotalDepth=refX+altX, minAltDepth=altX)
  
  # No mvr given
  if (!length(mvr)) return(mvr[0])

  submvr <- locateMTvariants(mvr, coding=coding)

  # If submvr is empty
  if (!length(submvr)) return(mvr[0])
  
  # mitochondrial genomic sequence
  # FIXME: may not want to do this
  # FIXME: else may want to filter

  data(rCRSeq, package="MTseeker")
  #altSeq <- DNAStringSet(replaceAt(rCRSeq[[1]], ranges(mvr), alt(mvr)))
  #names(altSeq) <- names(rCRSeq)
  #gr$varDNA <- getSeq(altSeq, gr)
  
  # use MT_CODE to translate results
  MT_CODE <- getGeneticCode("SGC1")
  gr$refSeq <- getSeq(rCRSeq, gr)
  
  # MT-RNR2 has an N in the sequence to allow for historic numbering to remain
  # So that has to be deleted before we can translate to get the AA seq
  if (genome(mvr) == "rCRS") {
    
    gr$refAA <- NA_character_
    
    # N is at bp 1737 of MT-RNR2
    gr["MT-RNR2"]$refSeq <- replaceAt(gr["MT-RNR2"]$refSeq, IRanges(1437,1437), "")
    gr$refAA <- suppressWarnings(translate(gr$refSeq, MT_CODE))
    
    # Now to keep numbering, we have to get the original reference sequence back
    gr["MT-RNR2"]$refSeq <- getSeq(rCRSeq, gr["MT-RNR2"])
    
  } else gr$refAA <- suppressWarnings(translate(gr$refSeq, MT_CODE))
  
  # This function make modifications based upon the reference
  # So for now the variant sequences are the same
  gr$varAA <- gr$refAA
  gr$varSeq <- gr$refSeq 
  
  gr$consequences <- NA_character_
  #gr$typeMut <- NA_character_
  gr$overlapGene <- NA_character_
  
  
  # Determine the type of mutation
  #submvr$typeMut <- NA_character_
  
  # Iterate through all of the affected genes
  for (i in 1:length(submvr)) {

    # Name of gene 
    g <- submvr[i]$gene
    
    # Range for the variant
    subir <- IRanges(submvr[i]$localStart, submvr[i]$localEnd)
    refSeq <- as.character(unlist(unname(extractAt(gr[g]$refSeq, subir))))
    
    # These sequences must be the same
    # Both are based off of the reference sequence
    if (refSeq != ref(submvr[i])) {
      if (g == "MT-ND6" && genome(submvr) == "rCRS") {
        gr[g]$consequences <- ""
        message("It is difficult to pinpoint an MT-ND6 variant start location, skipping . . .")
        next
      }
      else {
        message("Injecting variant incorrectly")
        show(submvr[i])
      }
      
    }
                        
    # Replace the part of the reference sequences with the variants
    # Then translates the variant sequences to get the amino acid sequence
    gr[g]$varSeq <- replaceAt(gr[g]$refSeq, subir, alt(submvr[i]))
    gr[g]$varAA <- suppressWarnings(translate(gr[g]$varSeq, MT_CODE))
    
    # AA sequence from the reference at the variant site
    # I was trying to be conservative on the codon position
    # Always assume the end of the AA sequence is end for getting the reference AA if that is the case
    endCodon <- submvr[i]$endCodon
    if (submvr[i]$endCodon > width(gr[g]$refAA)) endCodon <- width(gr[g]$refAA)
    orig <- extractAt(gr[g]$refAA, IRanges(submvr[i]$startCodon, endCodon))
    
    # If there is a deletion at the end of a gene
    # Must go here otherwise the endCodon value changes and you get the incorrect reference AA
    if (submvr[i]$endCodon > width(gr[g]$varAA)) {
      submvr[i]$endCodon <- width(gr[g]$varAA)
    }
    
    # Get the alternate AA sequences at the variant site
    altd <- extractAt(gr[g]$varAA, IRanges(submvr[i]$startCodon, submvr[i]$endCodon))
    gr[g]$consequences <- .flattenConsequences(orig, altd, submvr[i]$startCodon)
    
    gr[g]$overlapGene <- submvr[i]$overlapGene
    #if (refSeq != ref(submvr[i])) .compCheck(gr, g, submvr[i])
    # San check :)
    #if (!is.na(submvr[i]$overlapGene)) {
    #  show(names(submvr[i]))
    #  show(gr[g]$consequences)
    #  .compCheck(gr, g, submvr[i])
    #}
    
  } # for
  
  # later
  if (FALSE) {
    alignments <- pairwiseAlignment(gr$refAA, gr$varAA,
                                    substitutionMatrix = "BLOSUM50",
                                    gapOpening = 0, gapExtension = 8)
  }
  
  # Only return the genes that are affected by this variant
  return(gr[submvr$gene])
}

# helper function
.flattenConsequences <- function(orig, altd, startCodons) {
  
  .flat <- function(x) vapply(x[[1]], as.character, "N") 
  .prettify <- function(x) paste0(gsub(" ", "", x), collapse="")
  asdf <- DataFrame(refAA=.flat(orig), pos=startCodons, varAA=.flat(altd))
  csqs <- apply(subset(asdf, asdf$refAA != asdf$varAA), 1, .prettify)
  paste(csqs, collapse=",")
  
}

# helper fn
.compCheck <- function(gr, g, submvr) {
  
  print(as.character(names(submvr)))
  
  # DNA sequence ranges
  varRanges <- IRanges(submvr$localStart, width(gr[g]$varSeq))
  refRanges <- IRanges(submvr$localStart, width(gr[g]$refSeq))
  
  varDna <- unlist(unname(extractAt(gr[g]$varSeq, varRanges)))
  refDna <- unlist(unname(extractAt(gr[g]$refSeq, refRanges)))
  
  print("DNA sequences:")
  comp <- pairwiseAlignment(refDna, varDna)
  show(comp)
  
  # AA sequence ranges
  varRanges <- IRanges(submvr$startCodon, width(gr[g]$varAA))
  refRanges <- IRanges(submvr$startCodon, width(gr[g]$refAA))
  
  varAA <- unlist(unname(extractAt(gr[g]$varAA, varRanges)))
  refAA <- unlist(unname(extractAt(gr[g]$refAA, refRanges)))
  
  print("AA sequences")
  comp <- pairwiseAlignment(refAA, varAA)
  show(comp)
}

# May become obsolete with getProteinImpact function
.typeOfMutation <- function(submvr) {
  
  # Determine the overaching consequences (type of mutation)
  # SNP
  if (grepl(">", names(submvr[i]))) {
    if (gr[g]$consequences == "") submvr[i]$typeMut <- "synonymous"
    else if ("*" %in% unlist(altd)) submvr[i]$typeMut <- "nonsense"
    else submvr[i]$typeMut <- "missense"
  }
  
  # Insertions
  else if (grepl("ins", names(submvr[i]))) {
    if ( (nchar(alt(submvr[i])) - 1)  %% 3 == 0 ) {
      submvr[i]$typeMut <- "insertion"
      orig <- extractAt(gr[g]$refAA, IRanges(submvr[i]$startCodon, submvr[i]$startCodon))
    }
    else submvr[i]$typeMut <- "frameshift"
  }
  
  # Deletions
  else {
    if ( (nchar(ref(submvr[i])) - 1) %% 3 == 0) {
      submvr[i]$typeMut <- "deletion"
      altd <- extractAt(gr[g]$varAA, IRanges(submvr[i]$startCodon, submvr[i]$startCodon))
    } 
    else submvr[i]$typeMut <- "frameshift"
    
  }
  
  # If there is a frameshift mutation
  # Return the rest of the AA sequence as the consequence
  if (submvr[i]$typeMut == "frameshift") {
    
    # Want to return the rest of the sequence that has come super wonky
    orig <- extractAt(gr[g]$refAA, IRanges(submvr[i]$startCodon, width(gr[g]$refAA)))
    altd <- extractAt(gr[g]$varAA, IRanges(submvr[i]$startCodon, width(gr[g]$varAA)))
    gr[g]$consequences <- .flattenConsequences(orig, altd, submvr[i]$startCodon)
  }
  
  # Store the information
  gr[g]$typeMut <- submvr[i]$typeMut
  
}

.getAnno <- function(ref) {
  
  if (ref == "rCRS") {
    data("mtAnno.rCRS", package="MTseeker")
    anno <- mtAnno
  }
  
  return(anno)
}
