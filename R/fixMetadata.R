#' fix metadata for an MAlignmentsList or MVRangesList, if needed and possible.
#' 
#' This function is punishingly simple -- it just calls validMetadata() and
#' assigns attr('fixedMeta') from the result. If no fixing is required, the 
#' object's existing metadata is used. The object (with fixed metadata) is 
#' then returned. If automatic fixes are impossible, a message is generated.
#' 
#' @param x   an MAlignmentsList or MVRangesList
#'
#' @return    the object, with fixed metadata (if needed and possible) 
#' 
#' @export
fixMetadata <- function(x) {
  metadata(x) <- attr(validMetadata(x), "fixedMeta")
  return(x)
}
