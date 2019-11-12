#' convenience wrapper for dumping filtered variants directly to an indexed VCF
#'
#' @param vars 	an MVRanges or MVRangesList to filter, HGVS-annotate, and dump
#' @param file  the name of the resulting VCF file to produce (default is NULL)
#' @param ...   other arguments, passed to filterMTvars()
#'
#' @return	the filtered set of variants (as an MVRanges) that was exported 
#' 
#' @import 	VariantAnnotation
#'
#' @examples
#' \dontrun{
#' library(MTseekerData) 
#' filterMTvars(RONKSvariants$RO_1)
#' exportFilteredAsVcf(RONKSvariants$RO_1) # not written to a file 
#' }
#' @export 
exportFilteredAsVcf <- function(vars, file=NULL, ...) { 
  
  if (is(vars, "MVRangesList")) vars <- unlist(vars)
  res <- filterMTvars(vars, ...)
  if (!is.null(file)) writeVcf(res, file, index=TRUE)
  return(res)

}
