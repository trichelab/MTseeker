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
  gal <- GenomicAlignments:::GAlignmentsList(...)
  if (is.null(names(gal))) warning("This MAlignmentsList has no element names!")
  updateObject(new("MAlignmentsList", gal))
}


#' MAlignmentsList methods (centralized).
#'
#' @name      MAlignmentsList-methods
NULL


#' @rdname          MAlignmentsList-methods
#' 
#' @param object    an MAlignmentsList, usually lacking some cached information
#'
#' @return          an updated MAlignmentsList with cached metadata summaries
#'
#' @export
setMethod("updateObject", signature(object="MAlignmentsList"),
          function(object) { 
            if (!"Summary" %in% names(metadata(object))) {
              message("Caching element summaries in metadata(object)$Summary:") 
              metadata(object)$Summary <- Summary(object)
              message("Done.")
            }
            if (!"" %in% names(metadata(object))) {
              message("Caching BAM file listings in metadata(object)$BAMs:")
              metadata(object)$BAMs <- fileName(object)
              message("Done.")
            }
            return(object)
          }) 


#' @rdname          MAlignmentsList-methods
#'
#' @param x   an MAlignmentsList
#' 
#' @return    estimated coverage (numeric vector)
#'
#' @import    IRanges
#' 
#' @export
setMethod("coverage", signature(x="MAlignmentsList"),
          function(x) Summary(x)[, "genomeCoverage"])


#' @rdname          MAlignmentsList-methods
#'
#' @param x   an MAlignmentsList
#' 
#' @return    estimated coverage (numeric vector)
#'
#' @import    S4Vectors
#' 
#' @export
setMethod("runLength", signature(x="MAlignmentsList"),
          function(x) sapply(x, runLength))


#' @rdname          MAlignmentsList-methods
#' 
#' @param object  an MAlignmentsList
#' 
#' @return        BAM file summary for the MAlignmentsList 
#' 
#' @export
setMethod("fileName", signature(object="MAlignmentsList"),
          function(object) {
            if ("BAMs" %in% names(metadata(object))) {
              BAMs <- metadata(object)$BAMs
            } else {
              BAMs <- DataFrame(BAM=sapply(object, fileName),
                                readLength=unname(sapply(object, runLength)),
                                genome=unname(sapply(object, genome)))
              if (!is.null(names(object))) rownames(BAMs) <- names(object)
            } 
            if (!is.null(names(object))) {
              return(BAMs[names(object), ])
            } else { 
              return(BAMs)
            }
          })


#' @rdname          MAlignmentsList-methods
#' 
#' @param x    an MAlignmentsList
#' 
#' @return     a DataFrame
#'
#' @export
setMethod("Summary", signature(x="MAlignmentsList"),
          function(x) {
            if ("Summary" %in% names(metadata(x))) {
              dat <- metadata(x)$Summary
            } else {
              dat <- DataFrame(reads=sapply(x, length),
                               readLength=sapply(x, runLength),
                               genomeSize=sapply(x, runValue))
              if (!is.null(names(x))) rownames(dat) <- names(x)
              dat$genomeCoverage <- round(with(dat,
                                               (reads*readLength) / genomeSize))
            }
            # try to make the most of cached summaries... 
            if (!is.null(names(x))) dat <- dat[names(x), ] 
            return(dat)
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
            cat("Summary(object) returns a ")
            show(Summary(object))
            cat("-------\n", sep = "")
            cat("seqinfo: ", summary(seqinfo(object)), "\n", sep = "")
          })
