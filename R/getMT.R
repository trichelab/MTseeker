#' grab the mitochondrial reads from a BAM & estimate their fraction (of total)
#'
#' This purely a convenience function, and an incredibly convenient one at that.
#' FIXME: this function is still useful, but it desperately needs refactoring.
#' 
#' @param bam       a BAM filename, or DataFrame/SummarizedExperiment with $BAM
#' @param filter    filter on bam$mtCovg? (default is FALSE, don't filter)
#' @param plotMAPQ  plot distribution of mitochondrial mapping quality? (FALSE)
#' @param ...       additional args to pass scanBamParam(), such as mapqFilter
#'
#' @return          an MAlignments or MAlignmentsList object
#' 
#' @import GenomicAlignments
#' @import GenomeInfoDb
#' @import Rsamtools
#'
#' @importFrom Biostrings lcsuffix
#' @importFrom graphics par
#'
#' @examples
#' library(MTseekerData)
#' BAMdir <- system.file("extdata", "BAMs", package="MTseekerData")
#' patientBAMs <- paste0(BAMdir, "/", list.files(BAMdir, pattern="^pt.*bam$"))
#' (mal <- getMT(patientBAMs[1]))
#'
#' @export
getMT <- function(bam, filter=FALSE, plotMAPQ=FALSE, ...) {

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
      return(MAlignmentsList(lapply(bams, getMT)))
    } else { 
      message("No matching records.")
      return(NULL) 
    }
    # }}}
  }

  # for individual BAM files:
  bam <- as.character(bam) 
  bai <- paste0(bam, ".bai")
  if (!file.exists(bam)) stop(paste("Cannot find file", bam))
  if (!file.exists(bai)) indexBam(bam)
  bamfile <- BamFile(bam, index=bai, asMates=TRUE)
  chrMs <- grep("(chrM|MT)$", seqlevels(bamfile), value=TRUE)
  mtSeqLen <- seqlengths(bamfile)[chrMs] # GRCh37/38, hg38 & rCRS are identical

  # modified to handle xenograft BAMs if necessary 
  if (length(mtSeqLen) > 1) {
    message("This looks like a xenograft BAM...")
  } 
  mtGenome <- ifelse(any(mtSeqLen == 16569), "rCRS","other") # hg19 is "other"
  # Note: based on the BAM headers alone, we can't distinguish rCRS from RSRS
  if (mtGenome == "other") {
    message(bam, " may be aligned to hg19; length(", chrMs, ") is not 16569.")
    message("Patches for this and other human/nonhuman MT refs are welcome!")
    stop("Currently unsupported mitochondrial reference detected; exiting.")
  }

  idxStats <- idxstatsBam(bamfile)
  rownames(idxStats) <- idxStats$seqnames
  
  if (length(mtSeqLen) > 1) {
    message("Looks like a xenograft, splitting...")
    dropChars <- lcsuffix(chrMs[1], chrMs[2])
    genomes <- sub("_+", "", 
                   vapply(chrMs, 
                          function(x) substr(x, 1, nchar(x) - dropChars),
                          character(1)))
    idxStats$genome <- .getGenome(idxStats$seqnames, genomes)
    reads <- split(idxStats, idxStats$genome)
  } else { 
    idxStats$genome <- "rCRS"
    reads <- list(idxStats)
  }

  mtFrac <- vapply(reads, .getReadProportions, bam=bam, numeric(1)) 

  # this part can probably go away now that we pileup()
  mtRanges <- GRanges(chrMs, IRanges(1, mtSeqLen), strand="*")
  mtView <- BamViews(bam, bai, bamRanges=mtRanges)
  flags <- scanBamFlag(isUnmappedQuery=FALSE, 
                       isSecondaryAlignment=FALSE, 
                       isNotPassingQualityControls=FALSE, 
                       isDuplicate=FALSE) 
  mtParam <- ScanBamParam(flag=flags, what=c("seq","mapq"), which=mtRanges, ...)
  mtReads <- suppressWarnings(readGAlignments(mtView, 
                                              use.names=TRUE, # for revmapping
                                              param=mtParam)[[1]])
  
  if (length(mtFrac) > 1) { 
    seqlevels(mtReads) <- seqlevelsInUse(mtReads)
    mtReads <- split(mtReads, .getGenome(seqnames(mtReads)))
    for (i in names(mtReads)) {
      seqlevels(mtReads[[i]]) <- seqlevelsInUse(mtReads[[i]])
      mtChr <- grep("(MT|chrM)$", seqlevelsInUse(mtReads[[i]]), value=TRUE)
      isCircular(seqinfo(mtReads[[i]]))[mtChr] <- TRUE 
    }
    genome(mtReads) <- .getGenome(seqlevels(mtReads))
  } else {
    mtReads <- keepSeqlevels(mtReads, chrMs)
    isCircular(seqinfo(mtReads))[chrMs] <- TRUE 
    seqlevelsStyle(mtReads) <- "UCSC"
    genome(mtReads) <- mtGenome
  }
  attr(mtReads, "mtFrac") <- mtFrac

  if (plotMAPQ) { 
   if (length(mtFrac == 1)) {
     plot(density(mcols(mtReads)$mapq), type="h", col="red",
          xlab="MAPQ", ylab="fraction of reads with this MAPQ", 
          main=paste("Mitochondrial read mapping quality for\n", bam))
   } else { 
     par(mfrow=c(1, length(mtFrac)))
     for (i in names(mtReads)) {
       plot(density(mcols(mtReads[[i]])$mapq), type="h", col="red",
            xlab="MAPQ", ylab="fraction of reads with this MAPQ", 
            main=paste("MT read MAPQ for", seqlevelsInUse(mtReads[[i]])))
i    }
   }
  } 

  # MAlignments == wrapped GAlignments
  if (length(mtFrac) > 1) { 
    mal <- MAlignmentsList(lapply(mtReads, MAlignments, bam=bam))
  } else { 
    mal <- MAlignments(gal=mtReads, bam=bam)
  }
  attr(mal, "coverage") <- coverage(mal)
  return(mal)

}


# helper function:
.rCRSvsRSRS <- function(referenceSequence) { 
  switch(as.character(extractAt(referenceSequence, IRanges(523,524))[[1]]),
         "AC"="rCRS",
         "NN"="RSRS",
         "neither rCRS nor RSRS")
}


# helper function:
.getGenome <- function(x, genomes=c("GRCh38","GRCh37","hg38","GRCm38","mm10")) {
  if (length(x) > 1) vapply(x, .getGenome, genomes=genomes, character(1))
  else for (genome in genomes) if (grepl(genome, x)) return(genome) 
}


# helper function:
.getReadProportions <- function(idxStats, bam, chrM="(chrM|MT)") {

  MT <- grep(chrM, idxStats$seqnames, value=TRUE)
  mtReadCount <- idxStats[MT, "mapped"] 
  names(mtReadCount) <- MT
  nucReadCount <- sum(subset(idxStats, seqnames != MT)$mapped)
  mtFrac <- mtReadCount / (mtReadCount + nucReadCount)
  message(bam, " maps ", mtReadCount, " unique reads (~",
          round(mtFrac * 100, 1), "%) to ", MT, ".")
  return(mtFrac)

}
