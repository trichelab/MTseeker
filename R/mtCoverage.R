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
#' plotMtCoverage does what one might expect, and plots (read or call) coverage.
#'
#' plotStrandedMtCoverage does the same thing, but keeps track of which strand. 
#' 
#' @param x         an MAlignments or MVRanges
#' @param ...       other arguments to pass to GenomicAlignments::coverage()
#'
#' @aliases         plotMtCoverage
#' @aliases         plotStrandedMtCoverage
#'
#' @return          an RleList (or, invisibly for plot functions, a result list)
#' 
#' @import GenomicAlignments 
#'
#' @examples
#' library(MTseekerData)
#'
#' data(RONKSreads)
#' mtCoverage(RONKSreads$RO_1)
#' plotMtCoverage(RONKSreads$RO_1)
#'
#' data(RONKSvariants)
#' mtCoverage(RONKSvariants$RO_1)
#' plotMtCoverage(RONKSvariants$RO_1)
#'
#' par(mfrow=c(1,2))
#' plotMtCoverage(RONKSreads$NKS_1)
#' title("Read coverage for normal kidney sample 1") 
#' plotMtCoverage(RONKSreads$RO_1)
#' title("Read coverage for renal oncocytoma sample 1") 
#' 
#' par(mfrow=c(1,2))
#' plotStrandedMtCoverage(RONKSreads$NKS_1)
#' title("Stranded read coverage for normal kidney sample 1") 
#' plotStrandedMtCoverage(RONKSreads$RO_1)
#' title("Stranded read coverage for renal oncocytoma sample 1") 
#' @export
mtCoverage <- function(x, ...) { 
  if (is(x, "VRanges")) {
    res <- coverage(as(x, "VRanges"))
  } else if (is(x, "GAlignments")) {
    res <- coverage(as(x, "GAlignments"))
  } else { 
    stop("Don't know how to compute MT coverage for a ", class(x))
  }
  return(res)
} 


#' @rdname    mtCoverage
#' 
#' @import    circlize
#' 
#' @export
plotMtCoverage <- function(x, ...) { 

  CHR <- unique(seqnames(x)) 
  if (is(x, "MAlignments") | is(x, "MVRanges")) x <- mtCoverage(x)[[CHR]]
  message("Plotting mitochondrial coverage...")
  data(mtAnno.rCRS)
  anno <- mtAnno #.rCRS
  dat <- data.frame(name=names(anno), start=start(anno), end=end(anno))
  covg <- x
  ymax <- max(covg)
  circos.clear() 
  circos.par("clock.wise"=FALSE, "start.degree"=90, "gap.degree"=0, 
             "track.margin"=c(0.005, 0.005), "cell.padding"=c(0.005,0,0.005,0), 
             "points.overflow.warning"=FALSE)
  circos.genomicInitialize(data=dat, plotType=NULL, major.by=16569)

  # outside track: coverage on heavy strand
  colr <- colorRamp2(c(0,20,40,ymax), c("red","black","darkgreen","green"))
  p <- function(region, value, ...) { # {{{
    circos.genomicLines(region, value, col=colr(value), type="h")
  } # }}}
  circos.genomicTrack(.cov(covg, anno), track.height=0.15, 
                      ylim=c(0,ymax), bg.border=NA, panel.fun=p)
  circos.yaxis("right", labels.cex=0.5)

  # main track, gene names and such
  pfunGenes <- function(x, y) { # {{{
    xlim <- CELL_META$xlim
    ylim <- CELL_META$ylim
    gr <- anno[CELL_META$sector.index]
    ytop <- .height(gr) * ifelse(strand(gr) == "+", 1, 0)
    ybot <- .height(gr) * ifelse(strand(gr) == "-", -1, 0)
    lab <- ifelse(CELL_META$sector.index == "DLP", "CR", CELL_META$sector.index)
    circos.rect(xlim[1], ybot, xlim[2], ytop, col=gr$itemRgb)
    if (gr$region %in% c("rRNA", "coding", "D-loop") & gr$name != "HVR3") {
      circos.text(mean(xlim), .textloc(gr), lab, col="black", 
                  cex=.textcex(gr), font=.textbold(gr), facing="clockwise", 
                  niceFacing=TRUE)
    }
  } # }}}
  circos.track(panel.fun=pfunGenes, ylim=c(-1,1), track.height=0.5, 
               track.margin=c(0,0), bg.border=NA)

  res <- list(anno=dat, covg=x)
  invisible(res)
}


#' @rdname    mtCoverage
#' 
#' @import    circlize
#' 
#' @export
plotStrandedMtCoverage <- function(x, ...) { 

  CHR <- unique(seqnames(x)) 
  if (is(x, "MAlignments") | is(x, "MVRanges")) {
    x <- lapply(lapply(byStrand(x), mtCoverage), `[[`, CHR)
  }
  message("Plotting stranded mitochondrial coverage...")
  data(mtAnno.rCRS)
  anno <- mtAnno #.rCRS
  dat <- data.frame(name=names(anno), start=start(anno), end=end(anno))
  covgHeavy <- x[["+"]]
  covgLight <- x[["-"]]
  ymaxHeavy <- max(covgHeavy) 
  ymaxLight <- max(covgLight)
  ymax <- max(ymaxHeavy, ymaxLight)
  circos.clear() 
  circos.par("clock.wise"=FALSE, "start.degree"=90, "gap.degree"=0, 
             "track.margin"=c(0.005, 0.005), "cell.padding"=c(0.005,0,0.005,0), 
             "points.overflow.warning"=FALSE)
  circos.genomicInitialize(data=dat, plotType=NULL, major.by=16569)

  # outside track: coverage on heavy strand
  colr <- colorRamp2(c(0, 20, 40, ymax), c("red","black","darkgreen","green"))
  pH <- function(region, value, ...) { # {{{
    circos.genomicLines(region, value, col=colr(value), type="h")
  } # }}}
  circos.genomicTrack(.cov(covgHeavy, anno), track.height=0.15, 
                      ylim=c(0, ymax), bg.border=NA, panel.fun=pH)
  circos.yaxis("right", labels.cex=0.5)

  # main track, gene names and such
  pfunGenes <- function(x, y) { # {{{
    xlim <- CELL_META$xlim
    ylim <- CELL_META$ylim
    gr <- anno[CELL_META$sector.index]
    ytop <- .height(gr) * ifelse(strand(gr) == "+", 1, 0)
    ybot <- .height(gr) * ifelse(strand(gr) == "-", -1, 0)
    lab <- ifelse(CELL_META$sector.index == "DLP", "CR", CELL_META$sector.index)
    circos.rect(xlim[1], ybot, xlim[2], ytop, col=gr$itemRgb)
    if (gr$region %in% c("rRNA", "coding", "D-loop") & gr$name != "HVR3") {
      circos.text(mean(xlim), .textloc(gr), lab, col="black", 
                  cex=.textcex(gr), font=.textbold(gr), facing="clockwise", 
                  niceFacing=TRUE)
    }
  } # }}}
  circos.track(panel.fun=pfunGenes, ylim=c(-1,1), track.height=0.5, 
               track.margin=c(0,0), bg.border=NA)

  # inside track: coverage on light strand
  colr2 <- colorRamp2(c(0, -20, -40, -1*ymax), 
                      c("red","black","darkgreen","green"))
  pL <- function(region, value, ...) { # {{{
    circos.genomicLines(region, value, col=colr2(value),type="h",baseline="top")
  } # }}}
  circos.genomicTrack(.cov(covgLight, anno, direction="-"), track.height=0.15, 
                      ylim=c((-1*ymax), 0), bg.border=NA, panel.fun=pL)
  res <- list(anno=dat, covg=x)
  invisible(res)
}


# helper fn
.cov <- function(covg, anno, direction=c("+","-")) {
  direction <- match.arg(direction)
  gr <- GRanges(rep("chrM", length(runValue(covg))), 
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
