#' call mitochondrial variants from an MAlignments object 
#'
#' FIXME: figure out a way to reprocess extracted chrM/MT reads against rCRS,
#'        regardless of what reference they were originally aligned against.
#' 
#' @param mal         an MAlignments (or, potentially, an MAlignmentsList) 
#' @param ...         other optional arguments to pass to callVariants
#' @param parallel    try to call in parallel? (FALSE; buggy as hell at present)
#' @param verbose     be verbose? (FALSE; turn on for debugging purposes)
#'
#' @import gmapR
#' @import VariantTools
#' @import GenomicAlignments
#'
#' @export
callMT <- function(mal, ..., parallel=FALSE, verbose=FALSE){

  if (!is(mal, "MAlignments") & !is(mal, "MAlignmentsList")) {
    stop("callMT needs a MAlignments or MAlignmentsList to call variants.")
  } else if (is(mal, "MAlignmentsList")) { 
    message("Variant-calling an MAlignmentsList (this may melt your machine).")
    if (parallel == TRUE) { 
      stop("Parallel mtDNA variant calling is not working at the moment.") 
      return(MVRangesList(mclapply(mal, callMT)))
    } else { 
      return(MVRangesList(lapply(mal, callMT))) 
    }
  } else { 

    mtChr <- seqlevelsInUse(mal)
    mtGenome <- unique(genome(mal))
    gmapGenome <- paste("GmapGenome", "Hsapiens", mtGenome, sep=".")
    if (verbose) {
      message("Attempting to load ",gmapGenome," as reference genome index...") 
    }

    if (!requireNamespace(gmapGenome)) {
      if (mtGenome == "rCRS") {
        stop("MTseeker::indexMtGenome() will build a suitable ",gmapGenome,".")
      }
    }

    try(attachNamespace(gmapGenome), silent=TRUE)
    supportedSeqLevels <- seqlevels(get(gmapGenome)) 
    if (!mtChr %in% supportedSeqLevels) { 
      stop(mtChr, " is not among the supported seqlevels in ", gmapGenome, ".")
    }

    # for accounting purposes: 
    bamFile <- fileName(mal) 
    isCircular(mal) <- FALSE # for variant calling
    pars <- TallyVariantsParam(get(gmapGenome), 
                               minimum_mapq=20L,
                               high_base_quality=20L,
                               ignore_duplicates=TRUE, 
                               read_length=as(runLength(mal), "integer"),
                               which=as(seqinfo(mal)[mtChr], "GRanges"),
                               indels=TRUE)


    # needed for feedback during sometimes-slow bulk variant calling 
    message("Tallying variants for ", bamFile, "...") 
    tallied <- tallyVariants(bamFile, pars)

    if (verbose) message("QA'ing variants...") 
    QAed <- qaVariants(tallied)

    if (verbose) message("Calling ", bamFile, " vars against ", mtGenome, "...")
    filters <- VariantCallingFilters(...)
    res <- callVariants(QAed, calling.filters=filters)

    if (verbose) message("Filtering variants...")
    sampleNames(res) <- gsub(paste0(".", mtGenome), "", 
                             gsub("\\.bam", "", basename(bamFile)))
    res$PASS <- apply(softFilterMatrix(res), 1, all) == 1
    res$VAF <- altDepth(res) / totalDepth(res)
    genome(res) <- mtGenome

    if (verbose) message("Formatting variants...") 
    mvr <- MVRanges(res, coverage(mal))
    names(mvr) <- mtHGVS(mvr)
    return(mvr)
  }
}
