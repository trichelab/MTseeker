#' Import a VCF directly into MTseeker
#'
#' @name importVCF
#' 
#' @param x         character, list or dataframe for file name of VCF(and BAM)
#' @param gzip      Is the input gzipped? (FIXME: autodetect this)
#' @param biscuit   Is the input VCF from biscuit? (FALSE)
#' @param verbose   Print all messages during progress? (TRUE) 
#' 
#' @import VariantTools
#' @import VariantAnnotation
#' @importFrom utils read.table
#' 
#' @return an MVRangesList
#'
#' @examples 
#' 
#' library(MTseekerData)
#' VCFdir <- system.file("extdata", package="MTseekerData")
#' VCF <- file.path(VCFdir, list.files(VCFdir, pattern="*.vcf.gz$"))
#' mvr <- importVCF(VCF)[[1]]
#' mvr$FILTER <- ifelse(totalDepth(mvr) >= 20, "PASS", ".")
#' mvr$PASS <- mvr$FILTER == "PASS"
#' filterMT(mvr)
#' 
#' @export
importVCF <- function(x, gzip = TRUE, biscuit = FALSE, verbose = TRUE) {
  bam_files <- NULL
  if(is(x, 'character') || is(x, 'list')) {
    vcf_files <- x
  } else if(is(x, 'data.frame') || is(x, 'dataframe')) {
    if(!("VCF" %in% names(x))) {
      message("If given input is in some type of data.frame,")
      message("there must be column named \"VCF\" in it.")
      stop("Exiting.")
    } else {
      if("BAM" %in% names(x)) {
        if(verbose) message("Found column name \"BAM\" in input.")
        bam_files <- as.character(x$BAM)
      }
      vcf_files <- as.character(x$VCF)
    }
  } else {
    stop('Input type must be either character, list, or data.frame')
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
                       gzip = gzip,
                       biscuit = biscuit,
                       verbose = verbose)
    
    if(is.null(mvrl)) { mvrl <- t_mvrl }
    else { mvrl <- c(mvrl, t_mvrl)}
  }
  
  return(mvrl)
}

#' Convert a VCF to an MVRangesList object
#'
#' @name vcf2mvrl
#' 
#' @param vcf_filename   A filename/path to the corresponding VCF
#' @param bam_filename   A filename/path to the corresponding BAM
#' @param gzip           Is the input gzipped? (TRUE; really ought not matter)
#' @param biscuit        Is the input VCF from biscuit? (FALSE)
#' @param verbose        Print all messages during progress? (TRUE) 
#' 
#' @import VariantTools
#' @import VariantAnnotation
#' 
#' @return an MVRangesList
#' 
#' @examples 
#' 
#' library(MTseekerData)
#' VCFdir <- system.file("extdata", package="MTseekerData")
#' VCF <- file.path(VCFdir, list.files(VCFdir, pattern="*.vcf.gz$"))
#' mvr <- vcf2mvrl(VCF)[[1]]
#' mvr$FILTER <- ifelse(totalDepth(mvr) >= 20, "PASS", ".")
#' mvr$PASS <- mvr$FILTER == "PASS"
#' filterMT(mvr)
#' 
#' @export
vcf2mvrl <- function(vcf_filename, bam_filename = NULL, gzip = TRUE,
                     biscuit = FALSE, verbose = TRUE) {
  
  if(verbose) message("Loading data from ", vcf_filename)
  
  #check if the file is gzip'd
  if (gzip) {
    chr <- as.character(read.table(gzfile(vcf_filename), nrows = 1, sep=NULL, header = FALSE)$V1)
  } else {
    chr <- as.character(read.table(vcf_filename, nrows = 1, sep=NULL, header = FALSE)$V1)
    }
  
  if(startsWith(chr, 'chr')) {
    mt_chr <- 'chrM'
  } else {
    mt_chr <- 'MT'
  }
  
  rang <- GRanges(mt_chr, IRanges(1, 20000)) 
  ## Total length of chrM is smaller than 20K (16571-hg19 or 16569-rCRS)
  param <- ScanVcfParam(which=rang)
  
  mvrl <- NULL
  collVcfs <- suppressWarnings(readVcf(vcf_filename, param = param))
    
  if(verbose) message("Found ", ncol(collVcfs), " sample(s)")
  for(i in 1:ncol(collVcfs)) {
    collVcf <- collVcfs[, i]
    if (biscuit) {
      #filter out the CpG, CH, CHH, CHG methylation variant calls
      if (!("CX" %in% colnames(info(collVcf)))) stop("This doesn't look like a VCF from biscuit. Stopping.")
      collVcf <- collVcf[is.na(info(collVcf)$CX),]
    }
    if("DP" %in% names(geno(collVcf))) {
      dp <- unname(geno(collVcf)$DP[, 1])
      if ("AC" %in% names(geno(collVcf))) {
        alt_c <- geno(collVcf)$AC
        alt_c[is.na(alt_c == 0)] <- 0
        alt_c <- unlist(alt_c)
      } else { 
        alt_c <- unname(geno(collVcf)$AD[, , 2])
        ref_c <- dp - alt_c
      } 
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
                                   altDepth.field = "AC")
      vr$VAF <- alt_c / dp
    } else {
      message("Please look into your VCF file before running again")
      stop("Can't find [DP] column and/or [AD]/[AC] column...")
    }
 
        
    vr$PASS <- vr$FILTER == "PASS"
      
    if(ncol(collVcfs) == 1 && !is.null(bam_filename)) {
      cov <- .calc_cov(ranges_obj = vr, bam_file = bam_filename, chr = mt_chr, verbose = verbose)
    } else {
      cov <- 1
    }

    mvr <- keepSeqlevels(MVRanges(vr, coverage=cov), mt_chr, 
                         pruning.mode="coarse") 
    isCircular(mvr)[mt_chr]<- TRUE
    genome(mvr) <- "rCRS"
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

#helper function to calculate coverage
.calc_cov <- function(ranges_obj, bam_file, chr, verbose=TRUE) {
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
