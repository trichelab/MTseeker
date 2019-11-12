#' Mitochondrial genome coverage and plots for MAlignments or MVRanges objects
#' 
#' We co-opted the `coverage` method to retrieve approximate coverage depth
#' across the mitochondrial genome in MAlignments[List] and MVRanges[list],
#' so this function gives back what it was supposed to do (provide an Rle) 
#' and can allow for some subsetting (e.g. variant-supporting-read coverage)
#' that may be of interest when interpreting results. 
#'
#' The plotting functions can handle MAlignments or MVRanges objects directly.
#' 
#' plotMTCoverage does what one might expect, and plots (read or call) coverage.
#'
#' plotStrandedMTCoverage does the same thing, but keeps track of which strand. 
#' 
#' @param x         an MAlignments or MVRanges
#' @param ...       other arguments to pass to GenomicAlignments::coverage()
#'
#' @aliases         plotMTCoverage
#' @aliases         plotStrandedMTCoverage
#'
#' @return          an RleList (or, invisibly for plot functions, a result list)
#' 
#' @import GenomicAlignments 
#'
#' @examples
#' \dontrun{
#' library(MTseekerData)
#'
#' data(RONKSreads)
#' MTcoverage(RONKSreads$RO_1)
#' plotMTCoverage(RONKSreads$RO_1)
#'
#' data(RONKSvariants)
#' MTcoverage(RONKSvariants$RO_1)
#' plotMTCoverage(RONKSvariants$RO_1)
#'
#' par(mfrow=c(1,2))
#' plotMTCoverage(RONKSreads$NKS_1)
#' title("Read coverage for normal kidney sample 1") 
#' plotMTCoverage(RONKSreads$RO_1)
#' title("Read coverage for renal oncocytoma sample 1") 
#' 
#' par(mfrow=c(1,2))
#' plotStrandedMTCoverage(RONKSreads$NKS_1)
#' title("Stranded read coverage for normal kidney sample 1") 
#' plotStrandedMTCoverage(RONKSreads$RO_1)
#' title("Stranded read coverage for renal oncocytoma sample 1") 
#' }
#' @export
MTcoverage <- function(x, ...) { 
  if (is(x, "VRanges")) {
    res <- coverage(as(x, "VRanges"))
  } else if (is(x, "GAlignments")) {
    res <- coverage(as(x, "GAlignments"))
  } else { 
    stop("Don't know how to compute MT coverage for a ", class(x))
  }
  return(res)
} 


#' @rdname    MTcoverage
#' 
#' @import    circlize
#' 
#' @param x         an MAlignments or MVRanges
#' @param ref       string denoting reference sequence to use
#' 
#' @export
plotMTCoverage <- function(x, ref=c("rCRS", "NC_005089"), ...) { 

  if (length(ref) > 1) {
    message("You forgot to set a reference, current you can use: ", paste0(ref, sep=", "))
    stop()
  }
  
  if (is(x, "GAlignmentsList")) {
    
    if (length(x) > 1) {
      message("You can only plot the coverage of 1 sample at a time")
      stop()
    }
    
    else x <- unlist(x)
  }

  if (is(x, "GAlignments")) {
    bams <- unique(mcols(x)$bam)
    x <- MAlignments(x, bams)
  }
  
  circos.clear()
  
  CHR <- unique(seqnames(x)) 
  if (is(x, "MAlignments") | is(x, "MVRanges")) covg <- MTcoverage(x)[[CHR]]
  message("Plotting mitochondrial coverage...")
  
  ymax <- max(covg)
  
  # If you are plotting from the raw reads, the genome is not set yet
  genome(x) <- ref
  
  # Initialize the circos plot
  anno <- initMTcircos(x)

  # Plotting coverage part
  # outside track: coverage on heavy strand
  colr <- colorRamp2(c(0,20,40,ymax), c("red","black","darkgreen","green"))
  p <- function(region, value, ...) { # {{{
    circos.genomicLines(region, value, col=colr(value), type="h")
  } # }}}
  circos.genomicTrack(.cov(covg, anno), track.height=0.15, 
                      ylim=c(0,ymax), bg.border=NA, panel.fun=p)
  circos.yaxis("right", labels.cex=0.5)
  
  # Plot the genes of the circos plot
  genesMTcircos(x, anno, legends=T)
  
  # Store the information
  dat <- data.frame(name=names(anno), start=start(anno), end=end(anno))
  res <- list(anno=dat, covg=covg)
  invisible(res)
}


#' @rdname    MTcoverage
#' 
#' @import    circlize
#' 
#' @export
plotStrandedMTCoverage <- function(x, ...) { 
  
  if (genome(x) == "NC_005089") {
    stop("Unable to plot stranded coverage for mice right now")
  }
  
  circos.clear()
  
  anno <- initMTcircos(x)
  
  CHR <- unique(seqnames(x)) 
  if (is(x, "MAlignments") | is(x, "MVRanges")) {
    covg <- lapply(lapply(byStrand(x, anno), MTcoverage), `[[`, CHR)
  }
  message("Plotting stranded mitochondrial coverage...")
  
  covgHeavy <- covg[["heavy"]]
  covgLight <- covg[["light"]]
  ymaxHeavy <- max(covgHeavy) 
  ymaxLight <- max(covgLight)
  ymax <- max(ymaxHeavy, ymaxLight)

  # outside track: coverage on heavy strand
  colr <- colorRamp2(c(0, 20, 40, ymax), c("red","black","darkgreen","green"))
  pH <- function(region, value, ...) { # {{{
    circos.genomicLines(region, value, col=colr(value), type="h")
  } # }}}
  circos.genomicTrack(.cov(covgHeavy, anno), track.height=0.15, 
                      ylim=c(0, ymax), bg.border=NA, panel.fun=pH)
  circos.yaxis("right", labels.cex=0.5)

  # main track, gene names and such
  genesMTcircos(x, anno)

  # inside track: coverage on light strand
  colr2 <- colorRamp2(c(0, -20, -40, -1*ymax), 
                      c("red","black","darkgreen","green"))
  pL <- function(region, value, ...) { # {{{
    circos.genomicLines(region, value, col=colr2(value),type="h",baseline="top")
  } # }}}
  circos.genomicTrack(.cov(covgLight, anno, direction="-"), track.height=0.15, 
                      ylim=c((-1*ymax), 0), bg.border=NA, panel.fun=pL)
  
  dat <- data.frame(name=names(anno), start=start(anno), end=end(anno))
  res <- list(anno=dat, covg=covg)
  invisible(res)
}


# helper fn
.cov <- function(covg, anno, direction=c("+","-")) {
  direction <- match.arg(direction)
  gr <- GRanges(rep(seqnames(anno), length.out=length(runValue(covg)), 
                    length(runValue(covg))), 
                IRanges(start=start(covg), 
                        end=end(covg)),
                value=runValue(covg), 
                strand="*")
  if (direction == "-") gr$value <- -1 * gr$value
  ol <- findOverlaps(gr, anno) 
  gr$gene <- NA 
  gr[queryHits(ol)]$gene <- anno[subjectHits(ol)]$name
  as.data.frame(gr)[, c("gene","start","end","value")]
}
