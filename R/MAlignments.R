#' wraps a GAlignments with information about coverage and its target BAM file
#' 
#' The runLength slot stores readLength and named for historic reasons.
#' 
#' @import GenomicAlignments
#' @import Rsamtools
#' 
#' @exportClass MAlignments
setClass("MAlignments",
         representation(bam="character",
                        runLength="numeric"), 
         contains="GAlignments")


#' wrap a GAlignments for easier stats
#'
#' Normally the MAlignments constructor will be called by getMT(bam).
#' 
#' @rdname            MAlignments-methods
#' 
#' @param gal         a GAlignments
#' @param bam         a bam filename
#'
#' @return            an MAlignments 
#' 
#' @import            GenomicAlignments
#'
#' @examples
#' 
#' library(MTseekerData)
#' BAMdir <- system.file("extdata", "BAMs", package="MTseekerData")
#' BAMs <- paste0(BAMdir, "/", list.files(BAMdir, pattern=".bam$"))
#' mal <- getMT(BAMs[1])
#' class(mal) 
#' show(mal) 
#' 
#' @export
MAlignments <- function(gal, bam) { 
  if (!is(gal, "GAlignments")) stop("gal must be a GAlignments. Exiting.")
  if (length(bam) > 1) stop("bam must be a string naming a BAM file. Exiting.")
  new("MAlignments", gal, bam=bam, runLength=(median(qwidth(gal))-1))
}


#' MAlignments methods (centralized).
#'
#' Depending on how a generic was originally designated, the arguments to 
#' these methods can have various argument names, but all of them tend to 
#' take an MAlignments as their argument.
#'
#' @param x         an MAlignments
#' @param files     an MAlignments
#' @param object    an MAlignments
#' 
#' @return          various things, as appropriate to the methods
#'
#' @name            MAlignments-methods
NULL


#' @rdname          MAlignments-methods
#'
#' @export
setGeneric("genomeCoverage", 
           function(x) unname((length(x)*readLength(x))/genomeLength(x)))


#' @rdname    MAlignments-methods
#' 
#' @export
setMethod("Summary", signature(x="MAlignments"),
          function(x) {
            c(reads=length(x),
              readLength=readLength(x),
              genomeSize=genomeLength(x),
              genomeCoverage=genomeCoverage(x))
          })


#' @rdname    MAlignments-methods
#'
#' @export
setMethod("show", signature(object="MAlignments"),
          function(object) {
            callNextMethod()
            cat("  -------\n")
            cat(paste0("  ", round(genomeCoverage(object)), 
                       "x approximate read coverage."), "\n")
          })


#' @rdname    MAlignments-methods
#'
#' @export
setMethod("fileName", signature(object="MAlignments"),
          function(object) return(object@bam))


#' @rdname    MAlignments-methods
#'
#' @export
setMethod("scanBamHeader", signature(files="MAlignments"),
          function(files) {
            return(scanBamHeader(fileName(files)))
          })


#' @rdname    MAlignments-methods
#' 
#' @export 
setGeneric("readLength", 
           function(x) {
             if ("readLength" %in% slotNames(x)) {
               return(x@readLength)
             } else if ("runLength" %in% slotNames(x)) { 
               return(x@runLength)
             } else { 
               stop("Don't know how to get readLength for x.")
             }
           })


#' @rdname    MAlignments-methods 
#' 
#' @export
setMethod("readLength", "MAlignments", function(x) x@runLength)


#' @rdname    MAlignments-methods
#'
#' @export 
setGeneric("genomeLength", 
           function(x) sum(unname(seqlengths(x)[seqlevelsInUse(x)])))


#' @rdname    MAlignments-methods 
#' 
#' @export
setMethod("genomeLength", "MAlignments", 
           function(x) unname(seqlengths(x)[seqlevelsInUse(x)]))
