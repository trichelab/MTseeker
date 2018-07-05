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
    meta <- anno[CELL_META$sector.index]
    color <- meta$itemRgb
    xlim <- CELL_META$xlim
    ylim <- CELL_META$ylim
    ytop <- ifelse(strand(meta) == "+", ifelse(meta$region == "tRNA", .5, 1), 0)
    ybot <- ifelse(strand(meta) == "-", ifelse(meta$region == "tRNA",-.5,-1), 0)
    lab <- ifelse(CELL_META$sector.index == "DLP", "CR", CELL_META$sector.index)
    circos.rect(xlim[1], ybot, xlim[2], ytop, col=color)
    if (meta$region %in% c("rRNA", "coding", "D-loop")) {
      textloc <- ifelse(strand(meta) == "+", 0.5, -0.5)
      textcex <- ifelse(meta$name %in% c("MT-ATP8","MT-ATP6","HVR1","HVR2"), .7,
                        ifelse(meta$name %in% c("MT-ND3","MT-ND4L","MT-ND6",
                                                "MT-CO2","MT-CO3","MT-RNR1",
                                                "MT-RNR2"), .8, .9))
      circos.text(mean(xlim), textloc, sub("HVR", "HV", sub("HVR3", "", lab)), 
                  cex=textcex, col="black", facing="clockwise", niceFacing=TRUE)
    }
  }
  dat <- data.frame(name=names(anno), start=start(anno), end=end(anno))
  circos.par("clock.wise"=FALSE, start.degree=90, gap.degree=0)
  circos.genomicInitialize(data=dat, plotType=NULL, major.by=16569)
 
  # outside track: variant annotations (use granges() if variants is an MVRL)
  if (!is.null(variants)) message("Warning: variant support is very sucky")
  circos.track(track.height=0.1, ylim=c(0,1), bg.border=NA)
  
  # main track, gene names and such
  circos.track(panel.fun=pfun, ylim=c(-1,1), track.height=0.5, bg.border=NA)

  # inside track: plots of VAFs/indels and such? (also could use granges(MVRL))
  if (!is.null(matrices)) message("Warning: matrix support is very sucky too")
  circos.track(track.height=0.2, ylim=c(0,1), bg.border=NA)

  res <- list(anno=dat, pfun=pfun)
  invisible(res)
}
