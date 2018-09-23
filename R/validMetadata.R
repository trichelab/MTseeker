#' Ensure that the metadata caches in MAlignmentsLists and MVRangesLists are OK
#'
#' In order to avoid a lot of lengthy calculations, both MAlignmentsList and 
#' MVRangesList objects keep a cache of some relevant statistics and filenames
#' in their metadata slot. If these caches get stale, it can cause problems.
#' 
#' This function performs some sanity checks on the caches so that the above 
#' problems are unlikely to occur, provided that checkMetadataCache() is called
#' at sensible times. This function is NOT a replacement for validObject().
#' 
#' @param x           an MAlignmentsList or an MVRangesList 
#' 
#' @return            TRUE or FALSE (if FALSE, attr(res)$mismatches shows why)
#'
#' @export 
validMetadata <- function(x) { 
  if (!class(x) %in% c("MAlignmentsList", "MVRangesList")) {
    stop("MTseeker::validMetadata expects an MAlignmentsList or MVRangesList.")
  }
  msgs <- c() 
  if (!all(metadata(x)$summaryCols %in% colnames(metadata(x)$cache))) 
    msgs <- append(msgs, "Mandatory summary columns missing from cache.")
  if (!identical(rownames(metadata(x)$cache), names(x))) 
    msgs <- append(msgs, "Element names differ from cache rownames.")
  if (nrow(metadata(x)$cache) > length(x)) 
    msgs <- append(msgs, "Cache has more rows than the object has elements.")
  if (nrow(metadata(x)$cache) < length(x)) 
    msgs <- append(msgs, "Cache has fewer rows than the object has elements.")
  retval <- (length(msgs) == 0)
  if (length(msgs) > 0) {
    message("Errors found in cached metadata; attempting to fix on-the-fly.")
    attr(retval, "errors") <- msgs
  }
  if (!identical(rownames(metadata(x)$cache), names(x))) {
    nm0 <- setdiff(names(x), rownames(metadata(x)$cache))
    nm1 <- setdiff(rownames(metadata(x)$cache), names(x))
    attr(retval, "mismatches") <- DataFrame(
      name=c(nm0, nm1),
      inCache=c(rep(FALSE, length(nm0)), rep(TRUE, length(nm1))),
      inObject=c(rep(TRUE, length(nm0)), rep(FALSE, length(nm1)))
    )
  }
  fixedMeta <- metadata(x)
  if (length(msgs) > 0) {
    if (all(fixedMeta$summaryCols %in% colnames(fixedMeta$cache))) { 
      if (length(x) == nrow(metadata(x)$cache)) { 
        rownames(fixedMeta$cache) <- names(x)
      } else if (nrow(fixedMeta$cache) > length(x) &
                 all(names(x) %in% rownames(fixedMeta$cache))) {
        fixedMeta$cache <- fixedMeta$cache[names(x),]    
      }
    } else { 
      message("Object metadata has issues which cannot be fixed automatically.")
    }
  }
  attr(retval, "fixedMeta") <- fixedMeta
  return(retval)
}
