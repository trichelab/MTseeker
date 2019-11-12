#' process an MVRangesList for 20x depth and haplomask it
#' note that any elements with < 1 variant after filtering are tossed
#' 
#' @param   MVRL      An MVRangesList
#' @param   minDepth  minimum depth to retain (20) 
#' 
#' @import            VariantAnnotation
#' 
#' @export
depthFilterAndHaploMask <- function(MVRL, minDepth=20) {
  filtered <- MVRangesList(lapply(MVRL, subset, totalDepth >= minDepth))
  filtered <- filtered[which(sapply(filtered, length) > 0)] 
  haploMask(filtered)
}
