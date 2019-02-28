#' @import GenomicRanges
#' 
#' @param ranges_obj GRanges or MVRanges list or dataframe for file name of VCF(and BAM)
#' @param verbose    logical, if set to TRUE, print all messages during progressing.
#' @return           object as same as origina input
#'
#'
#' refered to http://haplogrep.uibk.ac.at/blog/rcrs-vs-rsrs-vs-hg19/
#' 
#' @export
lift_up_to_rCRS <- function(ranges_obj, verbose = T) {
  if(!(class(ranges_obj) == "GRanges" || class(ranges_obj) == "MVRanges")) {
    stop("Either GRanges or MVRanges must be provided as input")
  } 
  
  if(seqlengths(ranges_obj) == 16569) {
    stop("Input data already has 16569 length as same as chromosome M in version of rCRS")
  }
  
  if(verbose) message('Lifting up coordinates from HG19 base to rCRS base ... ')
  ir <- ranges(ranges_obj)
  
  ir1 <- ir[ir@start < 315]
  
  ir2 <- ir[ir@start >= 315 & (ir@start + ir@width - 1) < 3107]
  ir2 <- IRanges(start = (ir2@start - 2), end = (ir2@start + ir2@width - 1 - 2), names = names(ir2))
  
  ir3 <- ir[ir@start >= 3107 & (ir@start + ir@width - 1) < 16193]
  ir3 <- IRanges(start = (ir3@start - 1), end = (ir3@start + ir3@width - 1 - 1), names = names(ir3))
  
  ir4 <- ir[ir@start >= 16193]
  ir4 <- IRanges(start = (ir4@start - 2), end = (ir4@start + ir4@width - 1 - 2), names = names(ir4))
  
  ir <- c(ir1, ir2, ir3, ir4)
  
  res <- ranges_obj
  ranges(res) <- ir
  genome(res) <- "rCRS"
  seqlengths(res) <- 16569
  
  return(res)
}