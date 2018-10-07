#' mitochondrial genome coverage for an MAlignments or MVRanges object 
#' 
#' We co-opted the `coverage` method to retrieve approximate coverage depth
#' across the mitochondrial genome in MAlignments[List] and MVRanges[list],
#' so this function gives back what it was supposed to do (provide an Rle) 
#' and can allow for some subsetting (e.g. variant-supporting-read coverage)
#' that may be of interest when interpreting results. 
#' 
#' plotMtCoverage does what one might expect, and plots read coverage.
#' 
#' @param x         an MAlignments or MVRanges
#' @param ...       other arguments to pass to GenomicAlignments::coverage()
#'
#' @aliases         plotMtCoverage
#'
#' @return          an RleList
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
#' plotMtCoverage(RONKSreads$NKs_1)
#' title("Read coverage for normal kidney sample 1") 
#' plotMtCoverage(RONKSreads$RO_1)
#' title("Read coverage for renal oncocytoma sample 1") 
#' 
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
 
  if (is(x, "MAlignments") | is(x, "MVRanges")) x <- mtCoverage(x)
  message("Plotting mitochondrial coverage...")
  data(mtAnno.rCRS)
  anno <- mtAnno #.rCRS
  dat <- data.frame(name=names(anno), start=start(anno), end=end(anno))
  CHR <- grep("(chrM|MT)", names(x), value=TRUE) 
  covg <- x[[CHR]]
  ymax <- max(covg) 
  circos.clear() 
  circos.par("clock.wise"=FALSE, "start.degree"=90, "gap.degree"=0, 
             "track.margin"=c(0.005, 0.005), "cell.padding"=c(0.005,0,0.005,0), 
             "points.overflow.warning"=FALSE)
  circos.genomicInitialize(data=dat, plotType=NULL, major.by=16569)

  # outside track: coverage
  colr <- colorRamp2(c(0,20,40,max(covg)), c("red","black","darkgreen","green"))
  p <- function(region, value, ...) { # {{{
    circos.genomicLines(region, value, col=colr(value), type="h")
  } # }}}
  circos.genomicTrack(.cov(covg, anno), track.height=0.15, ylim=c(0,ymax),
                      bg.border=NA, panel.fun=p)
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

  res <- list(anno=dat)
  invisible(res)
}

# helper fn
.cov <- function(covg, anno) {
  gr <- GRanges(rep("chrM", length(runValue(covg))), 
                IRanges(start=start(covg), 
                        end=end(covg)),
                value=runValue(covg), 
                strand="*")
  gr$gene <- NA 
  ol <- findOverlaps(gr, anno) 
  gr[queryHits(ol)]$gene <- anno[subjectHits(ol)]$name
  as.data.frame(gr)[, c("gene","start","end","value")]
}
