#' simple helper function to split a *RangesList by the strand of its mt target
#' 
#' @param x     a *RangesList
#' @param anno  optional feature annotation, will use mtAnno.rCRS if NULL 
#'
#' @return      a list with elements of x overlapping features on each strand
#' 
#' @export 
byStrand <- function(x, anno=NULL) { 

  if (is.null(anno)) {
    data(mtAnno.rCRS)
    anno <- mtAnno
  }

  stranded <- split(anno, strand(anno))[c("+", "-")]
  names(stranded) <- c("heavy", "light")
  sapply(stranded, function(features) endoapply(x, subsetByOverlaps, features))

}
