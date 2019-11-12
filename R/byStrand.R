#' simple helper function to split a *RangesList by the strand of its mt target
#'
#' If presented with a GAlignments/MAlignments, this method will split the 
#' element by strand, i.e. + alignments and - alignments. Otherwise the method
#' retrieves ranges/variant calls that overlap genic elements on the heavy and
#' light strands of the mitochondrial genome. 
#' 
#' @param x     a *Ranges[List] or *Alignments[List]
#' @param anno  Annotation information for the reference genome
#'
#' @return      elements of x over features on each strand OR x split by strand 
#' 
#' @examples
#' \dontrun{
#' data(mtAnno.rCRS)
#' anno <- mtAnno
#' data(RONKSvariants, package="MTseekerData")
#' byStrand(RONKSvariants, anno)
#' }
#' @export 
byStrand <- function(x, anno) { 
  
  stranded <- split(anno, strand(anno))[c("+", "-")]
  names(stranded) <- c("heavy", "light")
  
  if (is(x, "MAlignments") | is(x, "GAlignments")) {
    split(x, strand(x))
  } else if (is(x, "GRanges") | is(x, "MVRanges")) {
    sapply(stranded, function(feats) subsetByOverlaps(x, feats))    
  } else { 
    sapply(stranded, function(feats) endoapply(x, subsetByOverlaps, feats))
  }

}
