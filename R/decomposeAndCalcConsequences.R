#' Decompose and annotate AA changes in MT variants
#'
#' @name decomposeAndCalcConsequences
#'
#' @param mvr    An MVRangesList or MVRanges object
#' @param AAchanges   Whether to annotate amino acid (AA) changes
#' @param parallel    Whether to run things in parallel
#' @param ...    Other arguments to pass to injectMTVariants
#'
#' @return    Annotated variants
#' 
#' @import GenomicRanges
#' @import parallel
#' @import VariantAnnotation
#' @import VariantTools
#'
#' @examples
#' 
#' library(MTseeker)
#' library(MTseekerData)
#' library(VariantTools)
#' 
#' #Set a really high depth filter
#' #This is just for an example and not something you'd use to filter real data
#' #Something like 10-20 reads is more reasonable
#' filters <- FilterRules(list(minTotalDepth = MinTotalDepthFilter(min.depth = 2000L)))
#' ronks_vars.anno <- RONKSvariants[1]
#' ronks_vars.anno <- MVRangesList(lapply(ronks_vars.anno, subsetByFilter, filters))
#' ronks_vars.anno <- decomposeAndCalcConsequences(ronks_vars.anno)
#' 
#' @export

decomposeAndCalcConsequences <- function(mvr, AAchanges=TRUE, parallel=FALSE, ...) {
  #this will decompose non-disjoint ranges for injectMTVariants()
  if (!class(mvr) %in% c("MVRanges", "MVRangesList")) stop("Input is not an MVRanges or MVRangesList.")
  #mvr.ovlps <- findOverlaps(mvr, type = "any")
  #get non-disjoint ranges
  #mvr.ovlps.nondisjoint <- mvr[queryHits(mvr.ovlps[queryHits(mvr.ovlps) != subjectHits(mvr.ovlps),]),]
  #keep disjoint ranges
  #mvr.ovlps.disjoint <- MVRanges(subsetByOverlaps(mvr, mvr.ovlps.nondisjoint, invert = TRUE))
  
  #run in parallel
  if (is(mvr, "MVRangesList") & parallel) {
    mvrl <- MVRangesList(mclapply(mvr, decomposeAndCalcConsequences))
    return(mvrl)
    }
  #run serially
  if (is(mvr, "MVRangesList")) {
    mvrl <- MVRangesList(lapply(mvr, decomposeAndCalcConsequences))
    return(mvrl)
    }

  #preprocess the variants
  mvr <- .getCoding(mvr, ...)
  
  if (length(mvr) == 0) {
    message("No variants overlapping coding space")
    return(mvr)
  }
  
  #add empty column for consequences
  mcols(mvr)$AAchange <- NA
  mcols(mvr)$impacted.gene <- NA
  
  if (isDisjoint(mvr)) {
    message("Found non-disjoint ranges in ", sampleNames(mvr)@values)
    message("Processing consequences...")
    if (AAchanges) {
      for (r in 1:length(mvr)) {
        con <- injectMTVariants(mvr[r])
        con.sub <- subset(con, mcols(con)$consequences != "")
        if (length(con.sub)) {
          if (length(mcols(con.sub)$consequences) > 1) {
            mcols(mvr)$AAchange[r] <- paste(mcols(con.sub)$consequences, collapse = ",")
          }
          if (length(mcols(con.sub)$synonym) > 1) {
            mcols(mvr)$impacted.gene[r] <- paste(mcols(con.sub)$synonym, collapse = ",")
          }
          else {
            mcols(mvr)$AAchange[r] <- mcols(con.sub)$consequences
            mcols(mvr)$impacted.gene[r] <- mcols(con.sub)$synonym 
          }
        } else {
          mcols(mvr)$impacted.gene[r] <- .getGeneImpacted(mvr[r])
          }
      }
    }
  } else {
    message("Processing consequences for ", sampleNames(mvr)@values)
    if (AAchanges) {
      for (r in 1:length(mvr)) {
        con <- injectMTVariants(mvr[r])
        con.sub <- subset(con, mcols(con)$consequences != "")
        if (length(con.sub)) {
          if (length(mcols(con.sub)$consequences) > 1) {
            mcols(mvr)$AAchange[r] <- paste(mcols(con.sub)$consequences, collapse = ",")
          }
          if (length(mcols(con.sub)$synonym) > 1) {
            mcols(mvr)$impacted.gene[r] <- paste(mcols(con.sub)$synonym, collapse = ",")
          }
          else {
            mcols(mvr)$AAchange[r] <- mcols(con.sub)$consequences
            mcols(mvr)$impacted.gene[r] <- mcols(con.sub)$synonym 
          }
        } else {
          mcols(mvr)$impacted.gene[r] <- .getGeneImpacted(mvr[r])
        }
      }
    }
  }
  
  return(mvr)
}

# helper function to subset ranges to just coding space
.getCoding <- function(mvr, gr=NULL, canon=.99, refX=1, altX=1) {
  # rCRS only, for the time being 
  stopifnot(unique(genome(mvr)) == "rCRS")
  
  # get mtGenes if needed 
  if (is.null(gr)) gr <- genes(mvr)
  stopifnot(unique(genome(gr)) == "rCRS")

  # subset the variants to those that overlap the target GRanges and are canon
  if (length(subsetByOverlaps(mvr, gr, type="within"))) {
    mvr <- subsetByOverlaps(mvr, gr, type="within")
    #drop anything that has an N base.. this also looks like a weird bug?
    mvr <- mvr[!grepl("N", mvr@alt),]
    #check again whether we've now cleared out all the variants
    #return an empty ranges if we have
    if (length(mvr) == 0) {
      mvr <- MVRanges(subsetByOverlaps(mvr, gr, type="within"))
    }
  } else { mvr <- MVRanges(subsetByOverlaps(mvr, gr, type="within")) }
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
