#' wraps a GAlignments with information about coverage and its target BAM file
#' 
#' @import GenomicAlignments
#' 
#' @exportClass MAlignments
setClass("MAlignments",
         representation(bam="character",
                        runLength="numeric"), 
         contains="GAlignments")


#' wrap a GAlignments for easier stats
#'
#' @param gal         a GAlignments
#' @param bam         a bam filename
#'
#' @return            an MAlignments 
#' 
#' @import            GenomicAlignments
#' 
#' @export
MAlignments <- function(gal, bam) { 
  if (!is(gal, "GAlignments")) stop("gal must be a GAlignments. Exiting.")
  if (length(bam) > 1) stop("bam must be a string naming a BAM file. Exiting.")
  new("MAlignments", gal, bam=bam, runLength=(median(qwidth(gal))-1))
}


#' MAlignments methods (centralized).
#'
#' @name      MAlignments-methods
NULL


#' @rdname    MAlignments-methods
#' 
#' @param x   an MAlignments
#' 
#' @return    estimated coverage (numeric)
#'
#' @export
setMethod("coverage", signature(x="MAlignments"),
          function(x) {
            unname( (length(x) * runLength(x)) / runValue(x) )
          })


#' @rdname    MAlignments-methods
#' 
#' @param x   an MAlignments
#' 
#' @return    read length
#'
#' @export
setMethod("runLength", signature(x="MAlignments"),
          function(x) {
            return(x@runLength)
          })


#' @rdname    MAlignments-methods
#' 
#' @param x   an MAlignments
#' 
#' @return    genome size 
#'
#' @export
setMethod("runValue", signature(x="MAlignments"),
          function(x) {
            unname(seqlengths(x)[seqlevelsInUse(x)])
          })


#' @rdname    MAlignments-methods
#' 
#' @param x   an MAlignments
#' 
#' @return    summary of x
#'
#' @export
setMethod("Summary", signature(x="MAlignments"),
          function(x) {
            c(reads=length(x),
              readLength=runLength(x),
              genomeSize=runValue(x),
              coverage=coverage(x))
          })


#' @rdname    MAlignments-methods
#'
#' @param x   an MAlignments
#' 
#' @export
setMethod("show", signature(object="MAlignments"),
          function(object) {
            callNextMethod()
            cat("  -------\n")
            cat(paste0("  ", round(coverage(object)), 
                       "x approximate read coverage."), "\n")
          })


#' @rdname    MAlignments-methods
#'
#' @param file      an MAlignments
#' 
#' @return          a BAM filename (string)
#'
#' @import          Rsamtools
#'
#' @export
setMethod("fileName", signature(object="MAlignments"),
          function(object) {
            return(object@bam)
          })


#' @rdname    MAlignments-methods
#'
#' @param bamRanges an MAlignments
#' 
#' @return          a BAM header
#'
#' @import          Rsamtools
#'
#' @export
setMethod("scanBamHeader", signature(files="MAlignments"),
          function(files) {
            return(scanBamHeader(fileName(files)))
          })
