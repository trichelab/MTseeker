#' wraps a GAlignmentsList (made up of MAlignments) for nicer viewing
#' 
#' This is largely deprecated with the move from gmapR to Rsamtools::pileup.
#' 
#' @import GenomicAlignments
#' @import S4Vectors
#' @import IRanges
#' 
#' @exportClass MAlignmentsList
setClass("MAlignmentsList", contains="GAlignmentsList")


#' wrap a GAlignmentsList for viewing
#'
#' Normally the MAlignmentsList constructor will be called by getMT. 
#' 
#' @rdname            MAlignmentsList-methods
#' 
#' @param ...         MAlignments
#'
#' @return            an MAlignments 
#' 
#' @import            GenomicAlignments
#' 
#' @examples
#'
#' library(MTseekerData)
#' BAMdir <- system.file("extdata", "BAMs", package="MTseekerData")
#' print(BAMdir)
#' BAMs <- paste0(BAMdir, "/", list.files(BAMdir, pattern=".bam$"))
#' print(BAMs)
#' mall <- getMT(BAMs[1]) # a xenograft
#' genomeLength(mall) 
#' readLength(mall) 
#' fileName(mall) 
#' scanBamHeader(mall) 
#'
#' @export
MAlignmentsList <- function(...) {

  gal <- GenomicAlignments::GAlignmentsList(...)
  if (is.null(names(gal))) warning("This MAlignmentsList has no element names!")
  new("MAlignmentsList", gal)

}


#' MAlignmentsList methods (centralized).
#'
#' Depending on how a generic was originally designated, the arguments to 
#' these methods can have various argument names, but all of them tend to 
#' take an MAlignmentsList as their argument.
#'
#' @param object    an MAlignmentsList
#' @param files     an MAlignmentsList
#' @param x         an MAlignmentsList
#' 
#' @return          various objects, as appropriate to the method 
#'
#' @name            MAlignmentsList-methods
NULL


#' @rdname          MAlignmentsList-methods
#' 
#' @export
setMethod("show", signature(object="MAlignmentsList"),
          function(object) {
            cat("MAlignmentsList object of length", length(object), "\n")
            cat("seqinfo: ", summary(seqinfo(object)), "\n", sep = "")
          })


#' @rdname          MAlignmentsList-methods
#'
#' @export
setMethod("fileName", "MAlignmentsList",
          function(object) vapply(object, fileName, character(1)))


#' @rdname          MAlignmentsList-methods
#'
#' @export
setMethod("scanBamHeader", "MAlignmentsList",
          function(files) vapply(files, scanBamHeader, list(1)))


#' @rdname          MAlignmentsList-methods 
#' 
#' @export
setMethod("readLength", "MAlignmentsList", 
          function(x) vapply(x, readLength, numeric(1)))


#' @rdname          MAlignmentsList-methods 
#' 
#' @export
setMethod("genomeLength", "MAlignmentsList", 
           function(x) vapply(x, genomeLength, numeric(1)))


#' @rdname          MAlignmentsList-methods
#'
#' @export
setMethod("genomeCoverage", "MAlignmentsList", 
          function(x) vapply(x, genomeCoverage, numeric(1)))
