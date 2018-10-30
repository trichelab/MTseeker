#' wraps a GAlignmentsList (made up of MAlignments) for nicer viewing
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
#' @rdname          MAlignmentsList-methods
#' 
#' @param ...         MAlignments
#'
#' @return            an MAlignments 
#' 
#' @import            GenomicAlignments
#' 
#' @examples
#' library(MTseekerData)
#' BAMdir <- system.file("extdata", "BAMs", package="MTseekerData")
#' BAMs <- paste0(BAMdir, "/", list.files(BAMdir, pattern=".bam$"))
#' targets <- data.frame(BAM=BAMs) 
#' rownames(targets) <- sapply(strsplit(basename(BAMs), "\\."), `[`, 1)
#' mall <- getMT(targets)
#' class(mall) 
#' show(mall) 
#'
#' @export
MAlignmentsList <- function(...) {

  # this must be done first: 
  mdat <- list()
  mdat$cache <- data.frame(BAM=sapply(..., fileName),
                           reads=sapply(..., length),
                           readLength=sapply(..., readLength), 
                           genomeSize=sapply(..., genomeLength), 
                           genome=unname(sapply(..., genome)),
                           nuclearReads=unname(sapply(..., attr, "nucReads")),
                           mitoVsNuclear=unname(sapply(..., attr, "mtVsNuc")))

  # options(stringsAsFactors) fix
  if (is.factor(mdat$cache$BAM)) {
    mdat$cache$BAM <- levels(mdat$cache$BAM)[mdat$cache$BAM] 
  }
  mdat$cache$genomeCoverage <- with(mdat$cache, 
                                    round((reads*readLength) / genomeSize))
  mdat$summaryCols <- c("reads", "readLength", 
                        "genomeSize", "genomeCoverage")

  # not relevant if only chrM reads  
  if (any(mdat$cache$nuclearReads > 0)) {
    mdat$summaryCols <- append(mdat$summaryCols, 
                               c("nuclearReads","mitoVsNuclear"))
  }

  # if cache is not prepped beforehand, this will clobber it: 
  gal <- GenomicAlignments::GAlignmentsList(...)
  # and no, I don't entirely understand why

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
#' Depending on how a generic was originally designated, the arguments to 
#' these methods can have various argument names, but all of them tend to 
#' take an MAlignmentsList as their argument.
#'
#' @param x         an MAlignmentsList
#' @param object    an MAlignmentsList
#' 
#' @return          various objects, as appropriate to the method 
#'
#' @name            MAlignmentsList-methods
NULL



#' @rdname          MAlignmentsList-methods
#'
#' @export
setMethod("genomeCoverage", signature(x="MAlignmentsList"),
          function(x) {
            covg <- metadata(x)$cache$genomeCoverage
            names(covg) <- names(x)
            return(covg)
          })


#' @export
setMethod("readLength", signature(x="MAlignmentsList"),
          function(x) {
            rl <- metadata(x)$cache$readLength
            names(rl) <- names(x)
            return(rl)
          })


#' @rdname          MAlignmentsList-methods
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
#' @export
setMethod("Summary", signature(x="MAlignmentsList"),
          function(x) {
            DataFrame(metadata(x)$cache[, metadata(x)$summaryCols])
          })


#' @rdname          MAlignmentsList-methods
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


