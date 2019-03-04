#' call mitochondrial variants against rCRS from an MAlignments[List] object 
#'
#' `callMTVars` is a helper function for callMT 
#'
#' FIXME: transition gmapR from import to suggestion
#' FIXME: use Rsamtools::pileup by default
#' FIXME: optional haplogroup masking?
#' 
#' @name  callMT
#'
#' @param mal         an MAlignments (or, potentially, an MAlignmentsList) 
#' @param ...         other arguments to pass to VariantTools::callVariants
#' @param parallel    try to run in parallel? (FALSE; this is super unstable)
#' @param cores       the number of cores to use for parallel processing
#' @param verbose     be verbose? (FALSE; turn on for debugging purposes)
#'
#' @return            an MVRanges (or, potentially, an MVRangesList) 
#' 
#' @import gmapR
#' @import VariantTools
#' @import GenomicAlignments
#'
#' @examples
#' 
#' library(MTseekerData)
#' BAMdir <- system.file("extdata", "BAMs", package="MTseekerData")
#' BAMs <- paste0(BAMdir, "/", list.files(BAMdir, pattern=".bam$"))
#' (mal <- getMT(BAMs[1]))
#' if (requireNamespace("GmapGenome.Hsapiens.rCRS", quietly=TRUE)) {
#'   (mvr <- callMT(mal))
#'   filt(snpCall(mvr))
#' } else { 
#'   message("You have not yet installed an rCRS reference genome.")
#'   message("Consider running the indexMTgenome() function to do so.")
#'   message("The RONKSvariants object in MTseekerData is a result of callMT.")
#' }
#' 
#' @export
callMT <- function (mal, ..., parallel = FALSE, cores = 1, verbose = FALSE) 
{
  if (!is(mal, "MAlignments") & !is(mal, "MAlignmentsList")) {
    stop("callMT needs an MAlignments or MAlignmentsList to call variants.")
  }
  else if (is(mal, "MAlignmentsList")) {
    if (unique(genome(mal)) != "rCRS") {
      stop("Genomes besides rCRS (aka GRCh37/GRCh38/hg38) are not supported.")
    }
    else {
      message("Variant-calling an MAlignmentsList (may melt your machine)...")
    }
    mtChr <- grep("(MT|chrM|NC_012920.1|rCRS)", seqlevelsInUse(mal), 
                  value = TRUE)
    
    #set the BPPARAM number of cores to 1 so parallel won't blow up
    #also has the side effect of making sure low depth samples also run
    BiocParallel::register(MulticoreParam(workers=1))
    
    if (parallel == TRUE) {
      warning("Parallel mtDNA variant calling is REALLY flaky at the moment.")
      options(mc.cores = cores)
      mvrl <- MVRangesList(mcmapply(callMTVars, BAM = metadata(mal)$cache$BAM, 
                                    SIZE = metadata(mal)$cache$readLength, GENOME = unique(genome(mal)), 
                                    CHR = rep(mtChr, length(mal))))
    }
    else {
      mvrl <- MVRangesList(mapply(callMTVars, BAM = metadata(mal)$cache$BAM, 
                                  SIZE = metadata(mal)$cache$readLength, GENOME = unique(genome(mal)), 
                                  CHR = rep(mtChr, length(mal))))
    }
    seqinfo(mvrl) <- seqinfo(mal)
    names(mvrl) <- names(mal)
    return(mvrl)
  }
  else if (is(mal, "MAlignments")) {
    mtChr <- grep("(MT|chrM|NC_012920.1|rCRS)", seqlevelsInUse(mal), 
                  value = TRUE)
    mvr <- callMTVars(BAM = fileName(mal), COV = genomeCoverage(mal), 
                      SIZE = readLength(mal), GENOME = unique(genome(mal)), 
                      CHR = mtChr)
    seqinfo(mvr) <- seqinfo(mal)
    return(mvr)
  }
}


#' @rdname callMT
#' 
#' @param BAM       the BAM filename (for callMTVars)
#' @param SIZE      the read length (for callMTVars; default is 75)
#' @param GENOME    the reference genome (for callMTVars; default is rCRS)
#' @param CHR       the mt contig name (for callMTVars; default is chrM)
#' @param COV       average read coverage (so we don't have to countBam)
#' 
#' @export
callMTVars <- function(BAM, SIZE=75, GENOME="rCRS", CHR="chrM", COV=NULL,
                       verbose=FALSE){

  gmapGenome <- paste("GmapGenome", "Hsapiens", GENOME, sep=".")
  if (verbose) {
    message("Attempting to load ", gmapGenome, " as reference genome index...") 
  }

  if (!requireNamespace(gmapGenome)) {
    if (GENOME == "rCRS") {
      stop("MTseeker::indexMTGenome() will build a suitable ",gmapGenome,".")
    }
  }

  # load the index (usually, rCRS skeleton key)
  try(attachNamespace(gmapGenome), silent=TRUE)
  # BE CAREFUL HERE: this can break variant calling 
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
  mvr <- keepSeqlevels(MVRanges(res, coverage=COV), CHR, pruning.mode="coarse") 
  isCircular(mvr)[CHR]<- TRUE
  names(mvr) <- MTHGVS(mvr)
  return(mvr)

}
