#' call mitochondrial variants from an MAlignments object 
#'
#' FIXME: figure out a way to reprocess extracted chrM/MT reads against rCRS,
#'        regardless of what reference they were originally aligned against.
#' 
#' @param mal         an MAlignments (or, potentially, an MAlignmentsList) 
#' @param ...         other optional arguments to pass to callVariants
#' @param rCRS        lift to rCRS if not hg38/GRCh38? (FALSE, and unsupported!)
#' @param verbose     be verbose? (FALSE; turn on for debugging purposes)
#'
#' @import gmapR
#' @import VariantTools
#' @import GenomicAlignments
#'
#' @export
callMT <- function(mal, ..., rCRS=FALSE, verbose=FALSE){

  if (rCRS == TRUE) stop("MTseeker does not currently support non-rCRS genomes")

  if (!is(mal, "MAlignments") & !is(mal, "MAlignmentsList")) {
    stop("callMT needs a MAlignments or MAlignmentsList to call variants.")
  } else if (is(mal, "MAlignmentsList")) { 
    message("Variant-calling an MAlignmentsList (this may melt your machine).")
    return(MVRangesList(lapply(mal, callMT)))
  }

  mtChr <- seqlevelsInUse(mal)
  mtGenome <- unique(genome(mal))
  gmapGenome <- paste("GmapGenome", "Hsapiens", mtGenome, mtChr, sep=".")
  if (verbose) {
    message("Attempting to load ", gmapGenome, " as reference genome index...") 
  }

  if (!requireNamespace(gmapGenome)) {
    if (mtGenome == "rCRS") {
      message("Check out ", system.file("extdata/rCRS.R", package="MTseeker"))
      message("if loading ", gmapGenome, " fails (you can build it yourself).")
    }
  }

  try(attachNamespace(gmapGenome), silent=TRUE)
  genome(mal) <- paste(mtGenome, mtChr, sep=".")
  isCircular(mal) <- FALSE # for variant calling
  pars <- TallyVariantsParam(get(gmapGenome), 
                             minimum_mapq=20L,
                             high_base_quality=20L,
                             ignore_duplicates=TRUE, 
                             read_length=as(runLength(mal), "integer"),
                             which=as(seqinfo(mal)[mtChr], "GRanges"),
                             indels=TRUE)
  
  if (verbose) message("Tallying variants...") 
  tallied <- tallyVariants(fileName(mal), pars)

  if (verbose) message("QA'ing variants...") 
  QAed <- qaVariants(tallied)

  if (verbose) message("Calling variants...") 
  filters <- VariantCallingFilters(...)
  res <- callVariants(QAed, calling.filters=filters)

  if (verbose) message("Filtering variants...")
  sampleNames(res) <- gsub(paste0(".", mtGenome), "", 
                           gsub("\\.bam", "", basename(fileName(mal))))
  res$PASS <- apply(softFilterMatrix(res), 1, all) == 1
  res$VAF <- altDepth(res) / totalDepth(res)
  genome(res) <- mtGenome

  if (verbose) message("Formatting variants...") 
  mvr <- MVRanges(res, coverage(mal))
  if (rCRS == TRUE) mvr <- rCRS(mvr)
  names(mvr) <- mtHGVS(mvr)
  return(mvr)
}
