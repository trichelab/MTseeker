#' plot a canonical human (or, in principle, any) mitochondrial genome 
#'
#' The default font sizes, orientations, etc. are optimized for a "cold" start;
#' if you want to fiddle with the details, crack open the code and modify it...
#' or alternatively, add sectors/dendrograms inside of this "framed" version.
#' 
#' @param anno      a GRanges (optional, defaults to mtAnno.rCRS if none given)
#' @param variants  optional MVRanges or MVRangesList to plot outside the circle
#' @param matrices  optional matrices of data to plot inside of the genes circle
#' 
#' @return          invisibly, a list: `anno` (data.frame) + `pfun` (panel.fun)
#'
#' @import circlize
#' 
#' @export 
mtCircos <- function(anno=NULL, variants=NULL, matrices=NULL) { 
  circos.clear() 
  data(mtAnno.rCRS)
  if (is.null(anno)) anno <- mtAnno #.rCRS
  pfun <- function(x, y) {
    xlim <- CELL_META$xlim
    ylim <- CELL_META$ylim
    gr <- anno[CELL_META$sector.index]
    ytop <- .height(gr) * ifelse(strand(gr) == "+", 1, 0)
    ybot <- .height(gr) * ifelse(strand(gr) == "-", -1, 0)
    lab <- ifelse(CELL_META$sector.index == "DLP", "CR", CELL_META$sector.index)
    circos.rect(xlim[1], ybot, xlim[2], ytop, col=gr$itemRgb)
    if (gr$region %in% c("rRNA", "coding", "D-loop") & gr$name != "HVR3") {
      circos.text(mean(xlim), .textloc(gr), lab, col="black", cex=.textcex(gr),
                  font=.textbold(gr), facing="clockwise", niceFacing=TRUE)
    }
  }
  dat <- data.frame(name=names(anno), start=start(anno), end=end(anno))
  circos.par("clock.wise"=FALSE, start.degree=90, gap.degree=0)
  circos.genomicInitialize(data=dat, plotType=NULL, major.by=16569)
 
  # outside track: variant annotations (use granges() if variants is an MVRL)
  if (!is.null(variants)) {
    gr2 <- granges(variants)
    coding <- predictCoding(variants)
    browser()
    circos.track(track.height=0.1, ylim=c(0,1), bg.border=NA)
  } else { 
    circos.track(track.height=0.1, ylim=c(0,1), bg.border=NA)
  }
  
  # main track, gene names and such
  circos.track(panel.fun=pfun, ylim=c(-1,1), track.height=0.5, bg.border=NA)

  # inside track: plots of VAFs/indels and such? (also could use granges(MVRL))
  if (!is.null(matrices)) message("Warning: matrix support is very sucky too")
  circos.track(track.height=0.2, ylim=c(0,1), bg.border=NA)

  res <- list(anno=dat, pfun=pfun)
  invisible(res)
}

# helper fn
.height <- function(gr) ifelse(gr$region %in% c("tRNA", "D-loop"), 0.5, 1)

# helper fn
.halfheight <- function(gr) gr$region %in% c("tRNA", "D-loop")

# helper fn
.textloc <- function(gr) {
  ifelse(.halfheight(gr), .25, .5) * ifelse(strand(gr) == "+", 1, -1)
}

# helper fn
.textbold <- function(gr) ifelse(gr$region %in% c("coding", "rRNA"), 3, 1)

# helper fn
.textcex <- function(gr) { 
  ifelse(gr$name %in% c("HVR1","HVR2"), .65,
         ifelse(gr$name %in% c("MT-ND3","MT-ND4L","MT-ND6","MT-CO2","MT-CO3",
                               "MT-RNR1","MT-RNR2","MT-ATP8","MT-ATP6"),.8,.9))
}
