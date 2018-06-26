#' wraps a GAlignmentsList (made up of MAlignments) for nicer viewing
#' 
#' @import GenomicAlignments
#' 
#' @exportClass MAlignmentsList
setClass("MAlignmentsList", contains="GAlignmentsList")


#' wrap a GAlignmentsList for viewing
#'
#' @param ...         MAlignments
#'
#' @return            an MAlignments 
#' 
#' @import            GenomicAlignments
#' 
#' @export
MAlignmentsList <- function(...) {

  # this must be done first: 
  mdat <- list()
  mdat$cache <- data.frame(BAM=sapply(..., fileName),
                           reads=sapply(..., length),
                           readLength=sapply(..., runLength), 
                           genomeSize=sapply(..., runValue), 
                           genome=unname(sapply(..., genome)))
  # options(stringsAsFactors) fix
  if (is.factor(mdat$cache$BAM)) {
    mdat$cache$BAM <- levels(mdat$cache$BAM)[mdat$cache$BAM] 
  }
  mdat$cache$genomeCoverage <- with(mdat$cache, 
                                    round((reads*readLength) / genomeSize))
  mdat$summaryCols <- c("reads", "readLength", "genomeSize", "genomeCoverage")

  # because otherwise this clobbers it: 
  gal <- GenomicAlignments:::GAlignmentsList(...)

  # name entries if possible
  if (is.null(names(gal))) {
    warning("This MAlignmentsList has no element names!")
  } else {
    rownames(mdat$cache) <- names(gal) 
  }

  # construct the object + its cache
  mal <- new("MAlignmentsList", gal)
  metadata(mal) <- mdat
  return(mal)

}


#' MAlignmentsList methods (centralized).
#'
#' @name            MAlignmentsList-methods
NULL


#' @rdname          MAlignmentsList-methods
#'
#' @param x         an MAlignmentsList
#' 
#' @return          estimated coverage (numeric vector)
#'
#' @import          IRanges
#' 
#' @export
setMethod("coverage", signature(x="MAlignmentsList"),
          function(x) {
            covg <- metadata(x)$cache$coverage
            names(covg) <- names(x)
            return(covg)
          })


#' @rdname          MAlignmentsList-methods
#'
#' @param x         an MAlignmentsList
#' 
#' @return          estimated coverage (numeric vector)
#'
#' @import          S4Vectors
#' 
#' @export
setMethod("runLength", signature(x="MAlignmentsList"),
          function(x) {
            rl <- metadata(x)$cache$readLength
            names(rl) <- names(x)
            return(rl)
          })


#' @rdname          MAlignmentsList-methods
#' 
#' @param object    an MAlignmentsList
#' 
#' @return          BAM file summary for the MAlignmentsList 
#' 
#' @export
setMethod("fileName", signature(object="MAlignmentsList"),
          function(object) {
            BAMs <- metadata(object)$cache$BAM
            names(BAMs) <- names(object)
            return(BAMs)
          }) 


#' @rdname          MAlignmentsList-methods
#' 
#' @param x         an MAlignmentsList
#' 
#' @return          a DataFrame
#'
#' @export
setMethod("Summary", signature(x="MAlignmentsList"),
          function(x) {
            metadata(x)$cache[, metadata(x)$summaryCols]
          })


#' @rdname          MAlignmentsList-methods
#'
#' @param objects   an MAlignmentsList
#' 
#' @export
setMethod("show", signature(object="MAlignmentsList"),
          function(object) {
            cat("MAlignmentsList object of length", length(object), "\n")
            cat("-------\n", sep = "")
            cat("Summary(object):\n")
            show(Summary(object))
            cat("-------\n", sep = "")
            cat("seqinfo: ", summary(seqinfo(object)), "\n", sep = "")
          })
