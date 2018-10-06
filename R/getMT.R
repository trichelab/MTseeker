#' grab the mitochondrial reads from a BAM and estimate their fraction
#'
#' nb. this could probably be done faster for a list of BAMs but it's not
#' nb. nb. this returns NuMt-depleted mitochondrial GenomicAlignments
#' nb. nb. nb. for the time being, this function ONLY supports rCRS/GRCh/hg38!
#' nb. nb. nb. nb. both chrM and mtGenome are now autodetected from BAM headers
#' nb. nb. nb. nb. nb. in the process of converting to rCRS, chrM becomes "chrM"
#' 
#' @param bam       a BAM filename, or DataFrame/SummarizedExperiment with $BAM
#' @param filter    filter on bam$mtCovg? (default is FALSE, don't filter)
#' @param parallel  load multiple BAMs in parallel, if possible? (FALSE)
#' @param plotMAPQ  plot distribution of mitochondrial mapping quality? (FALSE)
#' @param ...       additional args to pass scanBamParam(), such as mapqFilter
#'
#' @return          an MAlignments or MAlignmentsList object
#' 
#' @import GenomicAlignments
#' @import Rsamtools
#'
#' @examples
#' library(MTseekerData)
#' data(RONKSreads)
#' \donttest{
#'   RONKSreads <- getMT(BAMs)
#' } 
#' RONKSreads
#'
#' @export
getMT <- function(bam, filter=FALSE, parallel=FALSE, plotMAPQ=FALSE, ...) {

  # for lists/DataFrames/SEs of data:
  if (is(bam,"SummarizedExperiment")|is(bam,"DataFrame")|is(bam,"data.frame")) {
    # {{{
    if (is(bam, "DataFrame") | is(bam, "data.frame")) {
      if (! "BAM" %in% names(bam)) stop("data frame must have column `BAM`")
    } else {
      if (! "BAM" %in% names(colData(bam))) stop("SE must have colData()$BAM")
    }
    if (filter == TRUE) { 
      if (is(bam, "DataFrame") | is(bam, "data.frame")) {
        if (!"mtCovg" %in% names(bam)) stop("cannot filter on missing `mtCovg`")
      } else { 
        if (!"mtCovg" %in% names(colData(bam))) stop("missing colData()$mtCovg")
      }
      bam <- filterMT(bam) 
    }
    if (nrow(bam) > 0) {
      bams <- bam$BAM
      if (is(bam, "SummarizedExperiment")) {
        names(bams) <- colnames(bam)
      } else { 
        names(bams) <- rownames(bam)
      }
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
    # }}}
  }

  # for individual BAM files:
  bai <- paste0(bam, ".bai")
  if (!file.exists(bai)) indexBam(bam)
  bamfile <- BamFile(bam, index=bai, asMates=TRUE)
  chrM <- ifelse("chrM" %in% seqlevels(bamfile), "chrM", "MT")
  mtSeqLen <- seqlengths(bamfile)[chrM] # GRCh37/38, hg38 & rCRS are identical
  mtGenome <- ifelse(mtSeqLen == 16569, "rCRS", "other") # hg19 YRI is "other"
  # Note: based on the BAM headers (alone), we can't distinguish rCRS from RSRS
  if (mtGenome == "other") {
    message("MTseeker currently supports only rCRS-derived reference genomes.")
    message(bam, " may be aligned to hg19, as length(", chrM, ") is not 16569.")
    message("We find lifting hg19 (YRI) chrM to rCRS/RSRS can create problems.")
    message("(Patches for this and other human/nonhuman MT refs are welcome!)")
    stop("Currently unsupported mitochondrial reference detected; exiting.")
  }

  idxStats <- idxstatsBam(bamfile)
  rownames(idxStats) <- idxStats$seqnames
  mtReadCount <- idxStats[chrM, "mapped"] 
  nucReadCount <- sum(subset(idxStats, seqnames != chrM)$mapped)
  mtFrac <- mtReadCount / sum(idxStats[, "mapped"])
  message(bam, " maps ", mtReadCount, " unique reads (~",
          round(mtFrac * 100, 1), "%) to ", chrM, ".")

  mtRange <- GRanges(chrM, IRanges(1, mtSeqLen), strand="*")
  mtView <- BamViews(bam, bai, bamRanges=mtRange)
  flags <- scanBamFlag(isPaired=TRUE, 
                       isProperPair=TRUE, 
                       isUnmappedQuery=FALSE, 
                       hasUnmappedMate=FALSE, 
                       isSecondaryAlignment=FALSE, 
                       isNotPassingQualityControls=FALSE, 
                       isDuplicate=FALSE) 
  mtParam <- ScanBamParam(flag=flags, what=c("seq","mapq"), ...) 
  mtReads <- suppressWarnings(readGAlignments(mtView, param=mtParam)[[1]])
  attr(mtReads, "mtFrac") <- mtFrac
  mtReads <- keepSeqlevels(mtReads, chrM)
  isCircular(seqinfo(mtReads))[chrM] <- TRUE 
  seqlevelsStyle(mtReads) <- "UCSC"
  genome(mtReads) <- mtGenome

  if (plotMAPQ) {
    plot(density(mcols(mtReads)$mapq), type="h", col="red",
         xlab="MAPQ", ylab="fraction of reads with this MAPQ", 
         main=paste("Mitochondrial read mapping quality for\n", bam))
  }

  # MAlignments == wrapped GAlignments
  mal <- MAlignments(gal=mtReads, bam=bam)
  attr(mal, "coverage") <- coverage(mal)
  attr(mal, "nucReads") <- nucReadCount
  attr(mal, "mtReads") <- mtReadCount
  attr(mal, "mtVsNuc") <- mtReadCount/nucReadCount
  return(mal)

}


# helper function:
.rCRSvsRSRS <- function(referenceSequence) { 
  switch(as.character(extractAt(referenceSequence, IRanges(523,524))[[1]]),
         "AC"="rCRS",
         "NN"="RSRS",
         "neither rCRS nor RSRS")
}
