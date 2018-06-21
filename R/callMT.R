#' call mitochondrial variants against rCRS from an MAlignments[List] object 
#'
#' `callMtVars` is a helper function for callMT 
#'
#' FIXME: figure out a way to reprocess extracted chrM/MT reads against rCRS,
#'        regardless of what reference they were originally aligned against.
#' 
#' @name  callMT
#'
#' @param mal         an MAlignments (or, potentially, an MAlignmentsList) 
#' @param ...         other optional arguments to pass to callVariants
#' @param parallel    try to run in parallel? (FALSE; this is super unstable)
#' @param verbose     be verbose? (FALSE; turn on for debugging purposes)
#'
#' @import gmapR
#' @import VariantTools
#' @import GenomicAlignments
#'
#' @export
callMT <- function(mal, ..., parallel=FALSE, verbose=FALSE) {

  if (!is(mal, "MAlignments") & !is(mal, "MAlignmentsList")) {

    stop("callMT needs an MAlignments or MAlignmentsList to call variants.")

  } else if (is(mal, "MAlignmentsList")) { 
   
    if (unique(genome(mal)) != "rCRS") {
      stop("Genomes besides rCRS (aka GRCh37/GRCh38/hg38) are not supported.")
    } else { 
      message("Variant-calling an MAlignmentsList (may melt your machine)...")
    } 

    if (parallel == TRUE) { 
      warning("Parallel mtDNA variant calling is REALLY flaky at the moment.") 
      mvrl <- MVRangesList(mcmapply(callMtVars,
                                    BAM=metadata(mal)$cache$BAM,
                                    SIZE=metadata(mal)$cache$readLength,
                                    GENOME=metadata(mal)$cache$genome,
                                    CHR=rep(seqlevelsInUse(mal), length(mal))))
    } else { 
      mvrl <- MVRangesList(mapply(callMtVars,
                                  BAM=metadata(mal)$cache$BAM,
                                  SIZE=metadata(mal)$cache$readLength,
                                  GENOME=metadata(mal)$cache$genome,
                                  CHR=rep(seqlevelsInUse(mal), length(mal))))
    }
    seqinfo(mvrl) <- seqinfo(mal) 
    names(mvrl) <- names(mal)
    return(mvrl) 

  } else if (is(mal, "MAlignments")) { 

    mvr <- callMtVars(BAM=fileName(mal),
                      COV=coverage(mal),
                      SIZE=runLength(mal),
                      GENOME=unname(unique(genome(mal))),
                      CHR=seqlevelsInUse(mal))
    seqinfo(mvr) <- seqinfo(mal)
    return(mvr)

  }

}


#' @rdname callMT
#' 
#' @param BAM       the BAM filename (for callMtVars)
#' @param SIZE      the read length (for callMtVars; default is 75)
#' @param GENOME    the reference genome (for callMtVars; default is rCRS)
#' @param CHR       the mt contig name (for callMtVars; default is chrM)
#' @param COV       average read coverage (so we don't have to countBam)
#' 
#' @export
callMtVars <- function(BAM,SIZE=75,GENOME="rCRS",CHR="chrM",COV=NULL,verbose=F){

  gmapGenome <- paste("GmapGenome", "Hsapiens", GENOME, sep=".")
  if (verbose) {
    message("Attempting to load ", gmapGenome, " as reference genome index...") 
  }

  if (!requireNamespace(gmapGenome)) {
    if (GENOME == "rCRS") {
      stop("MTseeker::indexMtGenome() will build a suitable ",gmapGenome,".")
    }
  }

  # load the index (usually, rCRS skeleton key)
  try(attachNamespace(gmapGenome), silent=TRUE)
  whichRanges <- as(seqinfo(get(gmapGenome))[CHR], "GRanges")

  # for accounting purposes: 
  pars <- TallyVariantsParam(get(gmapGenome), 
                             minimum_mapq=20L,
                             high_base_quality=20L,
                             ignore_duplicates=TRUE, 
                             read_length=as(SIZE, "integer"),
                             which=whichRanges, 
                             indels=TRUE)

  # feedback during sometimes slow variant calling 
  message("Tallying variants for ", BAM, "...") 
  tallied <- tallyVariants(BAM, pars)
  if (verbose) message("QA'ing variants...") 
  QAed <- qaVariants(tallied)
  if (verbose) message("Calling ", BAM, " vars against ", GENOME, "...")
  filters <- VariantCallingFilters()
  res <- callVariants(QAed, calling.filters=filters)
  if (verbose) message("Filtering variants...")
  sampleNames(res) <- gsub(paste0(".", GENOME), "", 
                           gsub("\\.bam", "", basename(BAM)))
  res$PASS <- apply(softFilterMatrix(res), 1, all) == 1
  res$VAF <- altDepth(res) / totalDepth(res)
  genome(res) <- GENOME

  if (is.null(COV)) { 
    flags <- scanBamFlag(isPaired=TRUE, 
                         isProperPair=TRUE, 
                         isUnmappedQuery=FALSE, 
                         hasUnmappedMate=FALSE, 
                         isSecondaryAlignment=FALSE, 
                         isNotPassingQualityControls=FALSE, 
                         isDuplicate=FALSE) 
    mtParam <- ScanBamParam(flag=flags, 
                            mapqFilter=20, 
                            which=as(seqinfo(res)[CHR], "GRanges"),
                            what="qwidth") 
    bf <- BamFile(BAM, index=paste0(BAM, ".bai"))
    ct <- countBam(bf, param=mtParam)
    COV <- round((ct$records*SIZE)/ct$width)
  } 

  if (verbose) message("Formatting variants...") 
  mvr <- keepSeqlevels(MVRanges(res, COV), CHR, pruning.mode="coarse") 
  isCircular(mvr)[CHR]<- TRUE
  names(mvr) <- mtHGVS(mvr)
  return(mvr)

}
