#' Combine multiple MAlignments objects together
#'
#' @name combineMAlignments
#'
#' @param mall    An MAlignmentsList object
#' @param malToCombine    MAlignments object names to combine
#' @param replace    Whether to replace the merged object into the MAlignmentsList in place
#'
#' @return    Either a combined MAlignments object or an updated MAlignmentsList
#' 
#' @import GenomicAlignments
#'
#' @examples
#' 
#' library(MTseeker)
#' library(MTseekerData)
#' snames <- names(RONKSreads)[grepl("NK", names(RONKSreads))][1:3]
#' mergedMAlist <- combineMAlignments(RONKSreads, malToCombine=snames, replace=TRUE)
#' 
#' @export

combineMAlignments <- function(mall, malToCombine=NULL, replace=FALSE) {
  # check to make sure mal is an MAlignments object
  if (!is(mall, "MAlignmentsList")) {
    message("The input doesn't appear to be an MAlignmentsList object.")
    stop("You can generate one using getMT() or try and wrap it with MAlignmentsList().")
  }
  
  # check if any names were given
  if (is.null(malToCombine)) {
    message("Need to supply the names of MAlignments to combine.")
    stop("Can be accessed from the MAlignmentsList with names(yourMAlignmentsList).")
  }
  
  if (length(malToCombine) < 2) stop("Need to have more than one MAlignments object to combine.")
  
  # check if the user would like to modify the MAlignmentsList in place
  # FIXME: add in what is being replaced!
  if (replace) {
    message("Replacing the MAlignments objects with the combined object in place.")
  }
  
  # do the merging
  # this is not ideal...
  combinedMA <- mall[[malToCombine[1]]]
  
  # loop through and append
  for (i in 2:length(malToCombine)) {
    idx <- malToCombine[i]
    combinedMA <- append(combinedMA, mall[[idx]])
  }
  
  # replace the combined MAlignments object in place
  if (replace) {
    reducedMAlignments <- mall[malToCombine]
    reducedMAlignments <- validMetadata(reducedMAlignments)
    mall <- .replaceMAlignments(mall, malToCombine, combinedMA, reducedMAlignments)
  }
  
  if (replace) {
    return(mall)
  } else {
    return(combinedMA)
  }
}

# helper function
.replaceMAlignments <- function(MAlist, MAnames, combinedMA, reducedMAlignments) {
  
  # this approach is from rlist::list.remove()
  MAlist.names <- names(MAlist)
  m <- vapply(MAnames, "==", logical(length(MAlist)), MAlist.names)
  selector <- apply(m, 1L, any)
  MAlist <- MAlist[!selector]
  
  #update metadata
  message("Cleaning and recomputing the metadata...")
  combinedMdat <- .recomputeMAlignmentsMetadata(reducedMAlignments)
  rownames(combinedMdat$cache) <- paste0("combined_", MAnames[1])
  oldMdat <- metadata(validMetadata(MAlist))
  combinedMdat$cache <- rbind(oldMdat$cache,
                              combinedMdat$cache)
  #combine with old data
  MAlist <- lapply(MAlist, c)
  MAlist <- c(MAlist, list(combinedMA))
  names(MAlist)[length(MAlist)] <- paste0("combined_", MAnames[1])
  # this is the slow step...
  MAlist <- GenomicAlignments::GAlignmentsList(MAlist)
  MAlist <- new("MAlignmentsList", MAlist)
  metadata(MAlist) <- combinedMdat
  message("Done.")
  return(MAlist)
}

# helper function to update the metadata
.recomputeMAlignmentsMetadata <- function(mdat) {
  # recompute the metadata for the combined MAlignments object
  mdat.update <- list()
  mdat.update$cache <- data.frame(BAM = paste(metadata(mdat)$cache$BAM, collapse = ', '),
                            reads = sum(Summary(mdat)$reads),
                            readLength = mean(readLength(mdat)),
                            genomeSize = genomeLength(mdat),
                            genome = genome(mdat),
                            nuclearReads = sum(Summary(mdat)$nuclearReads),
                            mitoVsNuclear = sum(Summary(mdat)$reads) / sum(Summary(mdat)$nuclearReads),
                            genomeCoverage = sum(genomeCoverage(mdat)))
  mdat.update$summaryCols <- c("reads", "readLength", 
                               "genomeSize", "genomeCoverage",
                               "nuclearReads", "mitoVsNuclear")
  
  return(mdat.update)
}

# helper function to clean the cache/metadata
# .cleanMAlignmentsListCache <- function(MAlist) {
#   if (length(names(MAlist)) != nrow(metadata(MAlist)$cache)) {
#     # this is largely to prevent an error with names not existing
#     metadata(MAlist)$cache <- metadata(MAlist)$cache[rownames(metadata(MAlist)$cache) %in% names(MAlist),]
#   }
#   return(MAlist)
# }
