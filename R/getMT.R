#' grab the mitochondrial reads from a BAM and estimate their fraction
#'
#' nb. this could probably be done faster for a list of BAMs but it's not
#' nb. nb. this returns NuMt-depleted mitochondrial GenomicAlignments
#' FIXME: liftOver hg18/hg19-aligned results to rCRS for better variant calls
#' 
#' @param bam       a BAM filename, or a RangedSummarizedExperiment with $BAM
#' @param chrM      what the mitochondrial contig is called. Default is "chrM" 
#' @param mtGenome  what mitochondrial assembly was used (default is hg19) 
#' @param plotMAPQ  plot distribution of mitochondrial mapping quality? (FALSE)
#' @param filter    filter on colData(bam)$mtCovg? (default is TRUE, if an RSE)
#' @param parallel  load multiple BAMs in parallel, if possible? (FALSE)
#'
#' @import GenomicAlignments
#' @import Rsamtools
#'
#' @export
getMT <- function(bam, chrM="chrM", mtGenome="hg19", 
                  plotMAPQ=FALSE, filter=TRUE, parallel=FALSE){

  if (is(bam, "RangedSummarizedExperiment")) {
    if (! "BAM" %in% names(colData(bam))) stop("RSE must have colData()$BAM")
    if (filter == TRUE) { 
      if (!"mtCovg" %in% names(colData(bam))) stop("filter on colData()$mtCovg")
      bam <- filterMT(bam) 
    }
    if (nrow(bam) > 0) {
      bams <- bam$BAM
      names(bams) <- colnames(bam)
      # why lapply() by default? 
      # because most laptops will die otherwise!
      if (parallel == TRUE) {
        message("Loading multiple BAMs in parallel may kill your machine.")
        message("Set options('mc.cores') beforehand, and beware of swapping.")
        return(MAlignmentsList(mclapply(bams, getMT)))
      } else {
        return(MAlignmentsList(lapply(bams, getMT)))
      }
    } else { 
      message("No matching records.")
      return(NULL) 
    }
  }

  bai <- paste0(bam, ".bai")
  if (!file.exists(bai)) indexBam(bam)
  bamfile <- BamFile(bam, index=bai, asMates=TRUE)
  idxStats <- idxstatsBam(bamfile)
  rownames(idxStats) <- idxStats$seqnames
  mtFrac <- idxStats[chrM, "mapped"] / sum(idxStats[, "mapped"])
  message(bam, " has ~", round(mtFrac * 100), "% mitochondrial reads.")

  mtRange <- GRanges(chrM, IRanges(1, idxStats[chrM, "seqlength"]), "*")
  mtView <- BamViews(bam, bai, bamRanges=mtRange)
  if (!base::grepl(mtGenome, bam)) {
    message(mtGenome, " (supplied or default) isn't found in your bam filename")
  }
  flags <- scanBamFlag(isPaired=TRUE, 
                       isProperPair=TRUE, 
                       isUnmappedQuery=FALSE, 
                       hasUnmappedMate=FALSE, 
                       isSecondaryAlignment=FALSE, 
                       isNotPassingQualityControls=FALSE, 
                       isDuplicate=FALSE)

  mtParam <- ScanBamParam(flag=flags, what=c("seq","mapq")) 
  mtReads <- suppressWarnings(readGAlignments(mtView, param=mtParam)[[1]])
  attr(mtReads, "mtFrac") <- mtFrac
  genome(mtReads) <- mtGenome
  mtReads <- keepSeqlevels(mtReads, chrM)
  isCircular(seqinfo(mtReads)) <- TRUE 

  if (plotMAPQ) {
    plot(density(mcols(mtReads)$mapq), type="h", col="red",
         xlab="MAPQ", ylab="fraction of reads with this MAPQ", 
         main=paste("Mitochondrial read mapping quality for\n", bam))
  }

  # MAlignments == wrapped GAlignments
  mal <- MAlignments(gal=mtReads, bam=bam)
  attr(mal, "coverage") <- coverage(mal)
  return(mal)

}
