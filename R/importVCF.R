#' @import VariantTools
#' @import VariantAnnotation
#' 
#' @param x         character, list or dataframe for file name of VCF(and BAM)
#' @param verbose   logical, if set to TRUE, print all messages during progressing.
#' @return          an MVRangesList
#'
#' @export
importVCF <- function(x, verbose = T) {
  bam_files <- NULL
  if(is(x, 'character') || is(x, 'list')) {
    vcf_files <- x
  } else if(is(x, 'data.frame') || is(x, 'dataframe')) {
    if(!("VCF" %in% names(x))) {
      stop("If given input is in type of data.frame, there must be column name \"VCF\" in input.")
    } else {
      if("BAM" %in% names(x)) {
        if(verbose) message("Found column name \"BAM\" in input.")
        bam_files <- as.character(x$BAM)
      }
      vcf_files <- as.character(x$VCF)
    }
  } else {
    stop('Input must be given in type of either character, list or dataframe')
  }

  mvrl <- NULL
  bam_file <- NULL
  for(i in 1:length(vcf_files)) {
    vcf_file <- vcf_files[i]
    if(!file.exists(vcf_file)) { stop(paste0("Can't find ", vcf_file)) }
    if(!is.null(bam_files)) { 
      bam_file <- bam_files[i] 
      if(!file.exists(bam_file)) { stop(paste0("Can't find ", bam_file)) }
    }
    
    t_mvrl <- vcf2mvrl(vcf_filename = vcf_file,
                       bam_filename = bam_file,
                       verbose = verbose)
    
    if(is.null(mvrl)) { mvrl <- t_mvrl }
    else { mvrl <- c(mvrl, t_mvrl)}
  }
  
  return(mvrl)
}

vcf2mvrl <- function(vcf_filename, bam_filename = NULL, verbose = T) {
  if(verbose) message("Loading data from ", vcf_filename)
  
  chr <- as.character(read.table(vcf_filename, nrows = 1, sep=NULL, header = F)$V1)
  
  if(startsWith(chr, 'chr')) {
    mt_chr <- 'chrM'
  } else {
    mt_chr <- 'MT'
  }
  
  rang <- GRanges(mt_chr, IRanges(1, 20000)) ## Total length of chrM is smaller than 20K (16571-hg19 or 16569-rCRS)
  param <- ScanVcfParam(which=rang)
  
  mvrl <- NULL
  collVcfs <- suppressWarnings(readVcf(vcf_filename, param = param))
    
  if(verbose) message("Found ", ncol(collVcfs), " sample(s)")
  for(i in 1:ncol(collVcfs)) {
    collVcf <- collVcfs[, i]
      
    if(all(c("DP", "AC") %in% names(geno(collVcf)))){
      dp <- unname(geno(collVcf)$DP[, 1])
      alt_c <- geno(collVcf)$AC
      alt_c[is.na(alt_c == 0)] <- 0
      alt_c <- unlist(alt_c)
      ref_c <- dp - alt_c
        
      gr <- rowRanges(collVcf)
      gr$REF <- as.character(gr$REF)
      gr$ALT <- as.character(gr$ALT@unlistData)
      gr$DP <- dp
      gr$AC <- alt_c
      gr$RC <- ref_c
      gr$sampleNames <- header(collVcf)@samples[i]
        
      vr <- makeVRangesFromGRanges(gr, 
                                   ref.field = "REF", 
                                   alt.field = "ALT",
                                   totalDepth.field = "DP",
                                   refDepth.field = "RC",
                                   altDepth.field = "AC"
      )
      vr$VAF <- alt_c / dp
    } else {
      message("Please look into your VCF file before running again")
      stop("Can't find columns, either [DP] or [AC], or both ...")
    }
      
    vr$PASS <- vr$FILTER == "PASS"
      
    if(ncol(collVcfs) == 1 && !is.null(bam_filename)) {
      cov <- calc_cov(ranges_obj = vr, bam_file = bam_filename, chr = mt_chr, verbose = verbose)
    } else {
      cov <- 1
    }
      
    mvr <- keepSeqlevels(MVRanges(vr, coverage=cov), mt_chr, pruning.mode="coarse") 
    isCircular(mvr)[mt_chr]<- TRUE
    names(mvr) <- MTHGVS(mvr)
      
    if(seqlengths(mvr) == 16571) {
      mvr <- lift_up_to_rCRS(mvr, verbose = verbose)  
    }
    t_mvrl <- MVRangesList(mvr)
    names(t_mvrl) <- unique(gr$sampleNames)
      
    if(is.null(mvrl)) { mvrl <- t_mvrl }
    else { mvrl <- c(mvrl, t_mvrl) }
  }
    
  return(mvrl)
}

calc_cov <- function(ranges_obj, bam_file, chr, verbose=T) {
  if(verbose) message("Calculating coverage for region of chromosome M .. ")
  
  flags <- scanBamFlag(isPaired=TRUE, 
                       isProperPair=TRUE, 
                       isUnmappedQuery=FALSE, 
                       hasUnmappedMate=FALSE, 
                       isSecondaryAlignment=FALSE, 
                       isNotPassingQualityControls=FALSE, 
                       isDuplicate=FALSE) 
  mtParam <- ScanBamParam(flag=flags, 
                          mapqFilter=20, 
                          which=as(seqinfo(ranges_obj)[chr], "GRanges"),
                          what="qwidth") 
  bf <- BamFile(bam_file, index=paste0(bam_file, ".bai"))
  ct <- countBam(bf, param=mtParam)
  st <- scanBam(bf, param=mtParam)
  AVG_READ_SIZE <- mean(unname(st)[[1]]$qwidth)
  COV <- round((ct$records*AVG_READ_SIZE)/ct$width)
  
  message("Coverage: ", COV, "x")
  
  return(COV)
}
