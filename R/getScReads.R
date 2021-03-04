#' merge into MTseeker; could also do this with SmartSeq, etc.
#'
#' For 10X BAMs, cellTag is CB and fragTag is UB.
#' For Salmon BAMs, cellTag is CB and fragTag is UR.
#' For BioRad BAMs, cellTag is DB and fragTag is XB (I think).
#' Per-cell pileups can be generated directly from BAMs, though this is tricky.
#' I'm not 100 percent sure how to deal with UMI reconciliation above, though.
#'
#' @param   BAM         a BAM file or a BamFile object 
#' @param   cellTag     what tag to use for cell identification (CB)
#' @param   fragTag     what tag to use for fragment identification (UB)
#' @param   tagFilt     optional tag filter (eg. list(CB="ATCGATCGAT-1")) (NULL)
#' @param   dedupe      drop reads with missing or duplicate (CB,UB)? (FALSE)
#' @param   MTonly      only reaad in mitochondrial reads? (FALSE) 
#' @param   byCell      split by cell? (FALSE; only set TRUE for MT)
#' @param   byFrag      split by fragment? (FALSE; only set TRUE for MT)
#'
#' @return              a GAlignments or GAlignmentsList (if byCell/byFrag)
#'
#' @import GenomicAlignments
#' @import BiocParallel
#' @import Rsamtools
#' 
#' @export 
getScReads <- function(BAM, cellTag="CB", fragTag="UB", tagFilt=NULL, dedupe=FALSE, MTonly=FALSE, byCell=FALSE, byFrag=FALSE) {
  
  if (!is(BAM, "BamFile")) BAM <- BamFile(BAM)

  if (MTonly) {
    seqs <- grep("(chrM|MT)", seqlevels(BAM), value=TRUE) 
    if (length(seqs) < 1) stop("Error: no (chrM|MT) contigs found!")
    param <- ScanBamParam(which=as(seqinfo(BAM)[seqs], "GRanges"),
                          tagFilter=tagFilt, 
                          tag=c(cellTag, fragTag), 
                          what=c("seq"))
  } else { 
    param <- ScanBamParam(tagFilter=tagFilt,
                          tag=c(cellTag, fragTag), 
                          what=c("seq"))
  }

  message("Importing reads from ", BAM$path, "...", appendLF=FALSE) 
  gal <- readGAlignments(BAM, param=param) # single-end for now
  seqlevs <- seqlevelsInUse(gal)
  gal <- keepSeqlevels(gal, seqlevs) # massively simplifies scanning

  if (dedupe) gal <- .dedupe(gal, cellTag=cellTag, fragTag=fragTag)  
  
  message(length(gal), ifelse(dedupe, " deduplicated ", " "), 
          "fragments imported from ", basename(BAM$path), ".")  

  if (byCell & !byFrag) {
    .byCell(gal, cellTag=cellTag)
  } else if (byFrag & !byCell) {
    .byFrag(gal, fragTag=fragTag)
  } else if (byFrag & byCell) {
    .byCellFrag(gal, cellTag=cellTag, fragTag=fragTag)
  } else {
    return(gal)
  }
}
  

# helper fn
.dedupe <- function(gal, cellTag="CB", fragTag="UB") {
  message("Dropping duplicate (cell, fragment) pairs.")
  keep <- (!is.null(mcols(gal)[, cellTag]) & 
           !is.na(mcols(gal)[, cellTag]) & 
           !is.null(mcols(gal)[, fragTag]) & 
           !is.na(mcols(gal)[, fragTag]) & 
           !duplicated(paste0(mcols(gal)[, cellTag], "_", 
                              mcols(gal)[, fragTag])))
  return(gal[keep]) 
}


# helper fn 
.byCell <- function(gal, cellTag="CB") {
  gal <- keepSeqlevels(gal, seqlevelsInUse(gal))
  cells <- paste0(cellTag, ":", mcols(gal)[, cellTag])
  mcols(gal)[, cellTag] <- NULL 
  split(gal, cells)
}


# helper fn 
.byFrag <- function(gal, fragTag="UB") {
  gal <- keepSeqlevels(gal, seqlevelsInUse(gal))
  frags <- paste0(fragTag, ":", mcols(gal)[, fragTag])
  mcols(gal)[, fragTag] <- NULL 
  split(gal, frags)
}


# helper fn 
.byCellFrag <- function(gal, cellTag="CB", fragTag="UB") {
  gal <- keepSeqlevels(gal, seqlevelsInUse(gal))
  jointTag <- paste0(cellTag, ":", mcols(gal)[, cellTag], ";", 
                     fragTag, ":", mcols(gal)[, fragTag])
  mcols(gal)[, cellTag] <- NULL 
  mcols(gal)[, fragTag] <- NULL 
  split(gal, jointTag)
}


# helper fn
.fragTable <- function(gal, cellTag="CB", fragTag="UB", BPPARAM=SerialParam()) {
  gal <- keepSeqlevels(gal, seqlevelsInUse(gal))
  frags <- .fragsPerCell(gal, cellTag=cellTag, fragTag=fragTag, BPPARAM=BPPARAM)
  res <- data.frame(cell=names(frags), fragments=unlist(frags))
  res[order(res$fragments, decreasing=TRUE), ] 
}


# helper fn -- parallelize? 
.fragsPerCell <- function(gal, cellTag="CB", fragTag="UB", BPPARAM=SerialParam()) {
  gal <- keepSeqlevels(gal, seqlevelsInUse(gal))
  if (is(gal, "GAlignmentsList")) {
    bplapply(gal, .uniqueFrags, fragTag=fragTag, BPPARAM=BPPARAM)
  } else {
    bplapply(.byCell(gal, cellTag=cellTag), 
             .uniqueFrags, fragTag=fragTag, BPPARAM=BPPARAM)
  }
}


# helper fn
.uniqueFrags <- function(gal, fragTag="UB") {
  nlevels(factor(mcols(gal)[, fragTag]))
}


# helper fn
.getKnee <- function(fragTable) {
  cdf <- ecdf(fragTable$fragments)
  qs <- sapply(seq(0.01, 0.99, 0.01), quantile, x=cdf)
  cutoff <- names(which.max(diff(qs)/qs[-1]))
  round(qs[cutoff])
} 


# helper fn
.cellFilter <- function(gal, ft=NULL, cellTag="CB", fragTag="UB", 
                        BPPARAM=SerialParam()) {
  if (!is(gal, "GAlignmentsList")) {
    gal <- .byCell(gal, cellTag=cellTag)
  }
  if (is.null(ft)) {
    ft <- .fragsPerCell(gal, cellTag=cellTag, fragTag=fragTag, BPPARAM=BPPARAM)
  }
  gal[subset(ft$cell, ft$fragments > .getKnee(ft))]
}


# helper fn
.coverageByCell <- function(gal, BPPARAM=SerialParam()) { 
  gal <- keepSeqlevels(gal, seqlevelsInUse(gal))
  if (is(gal, "GAlignmentsList")) {
    res <- bplapply(gal, coverage, BPPARAM=BPPARAM)
  } else {
    res <- bplapply(.byCell(gal, cellTag=cellTag), coverage, BPPARAM=BPPARAM)
  }
  .txCoverage(res, BPPARAM=BPPARAM)
}


# helper fn 
.txCoverage <- function(covgs, BPPARAM=SerialParam()) { 
  unlist(bplapply(covgs, .weightedMeanCoverage, BPPARAM=BPPARAM))
}


# helper fn
.weightedMeanCoverage <- function(covg) {
  mean(sapply(covg, mean), weights=sapply(covg, length))
}
