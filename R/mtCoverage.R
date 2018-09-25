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
#' @export
mtCoverage <- function(x, ...) { 
  message("You may find plotMtCoverage useful to visualize mtCoverage results.")
  res <- switch(class(x), 
                MVRanges=coverage(as(x, "VRanges")),
                MAlignments=coverage(as(x, "GAlignments")))
  return(res)
} 


# helper fn
plotMtCoverage <- function(x, ...) { 
 
  if (is(x, "MAlignments")) x <- mtCoverage(x)

  data(mtAnno.rCRS)
  anno <- mtAnno #.rCRS
  CHR <- grep("(chrM|MT)", names(x), value=TRUE) 
  covg <- x[[CHR]]
  ymax <- max(covg) 

  pfun <- function(x, y) {
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
  }
  dat <- data.frame(name=names(anno), start=start(anno), end=end(anno))

  circos.clear() 
  circos.par("clock.wise"=FALSE, start.degree=90, gap.degree=0)
  circos.genomicInitialize(data=dat, plotType=NULL, major.by=16569)
  # outside track: coverage
  colr <- colorRamp2(c(0,20,40,max(covg)), c("red","black","darkgreen","green"))
  circos.genomicTrack(.cov(covg, anno), track.height=0.15, ylim=c(0, max(covg)),
                      bg.border=NA, panel.fun = function(region, value, ...) {
                        circos.genomicLines(region, value, col=colr(value), 
                                            type="h")
                      })
  circos.yaxis("right")
  # main track, gene names and such
  circos.track(panel.fun=pfun, ylim=c(-1,1), track.height=0.5, 
               track.margin=c(0,0), bg.border=NA)

  res <- list(anno=dat, pfun=pfun)
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
