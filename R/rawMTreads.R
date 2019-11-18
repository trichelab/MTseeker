#' Raw reads shown as GAlignments[List]
#'
#' @param bams    List of files
#'
#' @return        GAlignments or GAlignmentsList of raw reads
#'
#'
#' @examples 
#' 
#' library(MTseekerData)
#' BAMdir <- system.file("extdata", "BAMs", package="MTseekerData")
#' BAMs <- list.files(BAMdir, pattern="bam$")
#' BAMs <- paste0(BAMdir, "/", BAMs)
#' rawReads <- rawMTreads(BAMs)
#' rawReads_1 <- MAlignments(rawReads[[1]], BAMs[1])
#' plotMTCoverage(rawReads_1)
#' 
#' @export
rawMTreads <- function(bams) {
  
  # To return a list
  if (length(bams) > 1) {
    
    manyBams <- lapply(bams, rawMTreads)
    return(GAlignmentsList(manyBams))
  }
  
  else {
    sbp <- scanMT(bams)
    bamWhat(sbp) <- "seq"
    reads <- readGAlignments(file=bams, param=sbp)

    # Save the bam if they want to plot the coverage later
    mcols(reads)$bam <- bams
    return(reads)
  }
}
