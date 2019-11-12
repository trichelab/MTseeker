#' Figure out where the mitochondrial reads in a BAM are, so we can grab those.
#'
#' This purely a convenience function, and an incredibly convenient one at that.
#' 
#' @param bam       BAM (must be indexed!) file(s) or object(s) with a @bam slot
#' @param chrM      search pattern for mitochondrial contig ("(rCRS|chrM|MT)")
#' @param ...       additional args to pass scanBamParam(), such as mapqFilter
#'
#' @return          a ScanBamParam object for the BAM(s) to use with pileup() 
#' 
#' @import GenomicAlignments
#' @import GenomeInfoDb
#' @import Rsamtools
#' 
#' @examples
#' \dontrun{
#' library(MTseekerData)
#' BAMdir <- system.file("extdata", "BAMs", package="MTseekerData")
#' BAMs <- paste0(BAMdir, "/", list.files(BAMdir, pattern=".bam$"))
#'
#' sbp <- scanMT(BAMs[1])
#' show(sbp) 
#'
#' sbps <- scanMT(BAMs, mapqFilter=20)
#' show(sbps) 
#' }
#' @export
scanMT <- function(bam, chrM="(rCRS|chrM|MT)", ...) { 

  # for vectorized access, or just brute stupidity
  if (length(bam) > 1) return(sapply(bam, scanMT, ...)) 

  # for GAlignments/MAlignments objects (which will soon be phased out)
  if (is(bam, "GAlignments") & "bam" %in% slotNames(bam)) bam <- bam@bam
  bam <- as.character(bam) 
  if (!file.exists(bam)) stop(paste("Cannot find file", bam))
  bai <- paste(bam, "bai", sep=".")
  if (!file.exists(bai)) stop(paste("Cannot find index", bai))
  bh <- scanBamHeader(bam)[[1]]

  contigs <- names(bh$targets)
  mt <- grep(chrM, contigs, ignore.case=TRUE, value=TRUE)
  mtWhich <- GRanges(mt, IRanges(1, seqlengths(BamFile(bam))[mt]))
  mtFlags <- scanBamFlag(isUnmappedQuery=FALSE, hasUnmappedMate=FALSE, 
                         isSecondaryAlignment=FALSE, isDuplicate=FALSE, 
                         isNotPassingQualityControls=FALSE) 
  return(ScanBamParam(flag=mtFlags, which=mtWhich, what="strand", ...)) 

}
