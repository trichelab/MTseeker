#' call mitochondrial tumor/normal variants from an MAlignments[List] object 
#'
#' @name  callTNMT
#'
#' @param tmal        tumor MAlignments
#' @param nmal        normal MAlignments
#' @param ...         other optional arguments to pass to callVariants
#'
#' @import gmapR
#' @import VariantTools
#' @import GenomicAlignments
#'
#' @export
callTNMT <- function(tmal, nmal,...) { 

  if (!is(tmal, "MAlignments") | !is(nmal, "MAlignments")) { 
    stop("callTNMT needs a two MAlignments objects to call variants.")
  } 
  if (unique(genome(tmal)) != unique(genome(tmal))) {
    stop("callTNMT can only operate on BAMs with matching genomes (duh?)")
  }
  if (runLength(tmal) != runLength(nmal)) {
    stop("callTNMT is only meant for tumor-normal pairs with equal read length")
  }
  mtChr <- grep("(MT|chrM|NC_012920.1|rCRS)", seqlevelsInUse(tmal), value=TRUE)
  mvr <- callTnMtVars(BAMS=c(tumor=fileName(tmal), normal=filename(nmal)),
                      GENOME=unname(unique(genome(tmal))),
                      SIZE=runLength(tmal),
                      COV=coverage(tmal),
                      CHR=mtChr)
  seqinfo(mvr) <- seqinfo(tmal)
  return(mvr)

}


#' @rdname callTNMT
#' 
#' @param BAMS      the BAM filenames
#' @param SIZE      the read length
#' @param GENOME    the reference genome (default rCRS)
#' @param CHR       the mt contig name (default is chrM)
#' @param COV       the average read coverage
#' 
#' @export
callTnMtVars <- function(BAMS, SIZE, COV, GENOME="rCRS", CHR="chrM") {

  gmapGenome <- paste("GmapGenome", "Hsapiens", GENOME, sep=".")
  if (!requireNamespace(gmapGenome)) {
    if (GENOME == "rCRS") {
      stop("MTseeker::indexMtGenome() will build a suitable ",gmapGenome,".")
    }
  }
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
  message("Tallying variants for ", paste(BAMS, collapse=" vs. "), "...") 
  tallied <- callSampleSpecificVariants(BAMS["tumor"], BAMS["normal"], pars)
  QAed <- qaVariants(tallied)
  filters <- VariantCallingFilters()
  res <- callVariants(QAed, calling.filters=filters)
  sampleNames(res) <- gsub(paste0(".", GENOME), "", 
                           gsub("\\.bam", "", basename(BAMS["tumor"])))
  res$PASS <- apply(softFilterMatrix(res), 1, all) == 1
  res$VAF <- altDepth(res) / totalDepth(res)
  genome(res) <- GENOME
  mvr <- keepSeqlevels(MVRanges(res, coverage=COV), CHR, pruning.mode="coarse") 
  isCircular(mvr)[CHR]<- TRUE
  names(mvr) <- mtHGVS(mvr)
  return(mvr)

}
