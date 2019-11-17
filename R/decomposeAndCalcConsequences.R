#' Decompose and annotate AA changes in MT variants
#'
#' FIXME: this function is remarkably slow 
#' 
#' @name decomposeAndCalcConsequences
#'
#' @param mvr         An MVRangesList or MVRanges object
#' @param coding      Annotate only coding regions? (TRUE)
#' @param AAchanges   Annotate amino acid (AA) changes? (TRUE) 
#' @param ...         Other arguments to pass to injectMTVariants
#'
#' @return            Annotated variants, as an MVRanges
#' 
#' @import GenomicRanges
#' @import VariantAnnotation
#' @import VariantTools
#'
#' @examples
#' library(MTseeker)
#' library(MTseekerData)
#' library(VariantTools)
#' 
#' # Set a really high depth filter
#' # Just an example, not something you'd use to filter real data
#' # Something like 20 reads is more reasonable (gives Sanger-like error rates)
#' filters <- 
#'   FilterRules(list(minTotalDepth=MinTotalDepthFilter(min.depth=2e3L)))
#' ronks_vars.anno <- subsetByFilter(RONKSvariants[[1]], filters)
#' decomposeAndCalcConsequences(ronks_vars.anno)
#' 
#' @export
decomposeAndCalcConsequences <- function(mvr, coding=TRUE, AAchanges=TRUE, ...){
  
  #this will decompose non-disjoint ranges for injectMTVariants()
  if (!class(mvr) %in% c("MVRanges", "MVRangesList")) stop("Input is not an MVRanges or MVRangesList.")
  
  if (is(mvr, "MVRangesList")) {
    MVRangesList(lapply(mvr, decomposeAndCalcConsequences, coding=coding, ...))
  }

  # Save the coverage information for later when we combine mvr 
  # and overlapMvr at the end
  covg <- genomeCoverage(mvr)

  #preprocess the variants
  if (coding) mvr <- .getCoding(mvr, ...)

  if (length(mvr) == 0) {
    message("No variants found within the coding region, returning empty MVRanges")
    return(mvr)
  }
  
  #add empty column for consequences
  mcols(mvr)$AAchange <- NA
  #mcols(mvr)$typeMut <- NA
  
  mcols(mvr)$impacted.gene <- NA
  mcols(mvr)$overlapGene <- NA
  
  # Store overlapping variants here
  # Even though there will be doubles of variants
  # At least you can see the potential depending on which gene it effects
  overlapMvr <- mvr[0]
  
  if (!isDisjoint(mvr)) {
    
    #message("Found non-disjoint ranges in ", sampleNames(mvr)@values)
    #message("Processing consequences...")
    
    if (AAchanges) {

      for (r in 1:length(mvr)) {
        
        con <- injectMTVariants(mvr[r], coding=coding, ...)
        
        if (length(con) == 2) {
          
          newMvr <- mvr[r]
          
          mcols(newMvr)$AAchange <- mcols(con)$consequences[2]
          if (mcols(newMvr)$AAchange == "") mcols(newMvr)$AAchange <- NA_character_
          #mcols(newMvr)$typeMut <- mcols(con)$typeMut[2]
          
          mcols(newMvr)$impacted.gene <- mcols(con)$synonym[2]
          mcols(newMvr)$overlapGene <- mcols(con)$overlapGene[2]
          
          overlapMvr <- append(overlapMvr, newMvr)
          
        } else {
          
          mcols(mvr)$AAchange[r] <- mcols(con)$consequences
          if (mcols(mvr)$AAchange[r] == "") mcols(mvr)$AAchange[r] <- NA_character_
          #mcols(mvr)$typeMut[r] <- mcols(con)$typeMut
          
          mcols(mvr)$impacted.gene[r] <- mcols(con)$synonym
          mcols(mvr)$overlapGene[r] <- mcols(con)$overlapGene
          
        }
      }
      
    } # AAchanges
    
  }
  
  else {
    #message("Processing consequences for ", sampleNames(mvr)@values)
    if (AAchanges) {
      
      for (r in 1:length(mvr)) {
        con <- injectMTVariants(mvr[r], coding=coding, ...)
        
        if (length(con) == 2) {

          newMvr <- mvr[r]
          
          mcols(newMvr)$AAchange <- mcols(con)$consequences[2]
          #mcols(newMvr)$typeMut <- mcols(con)$typeMut[2]
          
          mcols(newMvr)$impacted.gene <- mcols(con)$synonym[2]
          mcols(newMvr)$overlapGene <- mcols(con)$overlapGene[2]
          
          overlapMvr <- append(overlapMvr, newMvr)
          
        } else {
          
          mcols(mvr)$AAchange[r] <- mcols(con)$consequences
          #mcols(mvr)$typeMut[r] <- mcols(con)$typeMut
          
          mcols(mvr)$impacted.gene[r] <- mcols(con)$synonym
          mcols(mvr)$overlapGene[r] <- mcols(con)$overlapGene
          
        }
      }
      
      
    } ## AAchanges
  }
  

  mvr <- sort(MVRanges(c(mvr, overlapMvr), coverage = covg), ignore.strand=T)
  return(mvr)
}

# helper function to subset ranges to just coding space
.getCoding <- function(mvr, gr=NULL, canon=.99, refX=1, altX=1) {

  # rCRS only, for the time being 
  stopifnot(unique(genome(mvr)) == "rCRS")
  
  # get mtGenes if needed 
  if (is.null(gr)) gr <- genes(mvr)
  stopifnot(unique(genome(gr)) == "rCRS")

  mvr <- MVRanges(subsetByOverlaps(mvr, gr, type="within"), coverage = genomeCoverage(mvr))
  
  # subset the variants to those that overlap the target GRanges and are canon
  if (length(mvr)) {
    
    #drop anything that has an N base.. this also looks like a weird bug?
    mvr <- mvr[!grepl("N", mvr@alt),]
    
    #check again whether we've now cleared out all the variants
    #return an empty ranges if we have
    if (length(mvr) == 0) {
      mvr <- MVRanges(subsetByOverlaps(mvr, gr, type="within"))
    }
  } 
  return(mvr)
}

# helper function to pull of which gene the variant is impacting
.getGeneImpacted <- function(mvr, gr=NULL) {
  # rCRS only, for the time being 
  stopifnot(unique(genome(mvr)) == "rCRS")
  
  # get mtGenes if needed 
  if (is.null(gr)) gr <- genes(mvr)
  stopifnot(unique(genome(gr)) == "rCRS")
  
  gene.name <- mcols(subsetByOverlaps(gr, mvr))$synonym
  if (length(gene.name) > 1) {
    gene.name <- paste(gene.name, collapse = ",")
  }
  
  return(gene.name)
}
