#' plot a canonical human (or, in principle, any) mitochondrial genome 
#'
#' The default font sizes, orientations, etc. are optimized for a "cold" start;
#' if you want to fiddle with the details, crack open the code and modify it...
#' or alternatively, add sectors/dendrograms inside of this "framed" version.
#' 
#' @param anno      a GRanges (optional, defaults to mtAnno.rCRS if none given)
#' @param variants  an MVRanges or MVRangesList (optional, and poorly supported)
#' 
#' @return          invisibly, a list: `anno` (data.frame) + `pfun` (panel.fun)
#'
#' @import circlize
#' 
#' @export 
mtCircos <- function(anno=NULL, variants=NULL) { 
  circos.clear() 
  if (is.null(anno)) {
    data(mtAnno.rCRS)
    anno <- mtAnno #.rCRS
  }
  pfun <- function(x, y) {
    meta <- anno[CELL_META$sector.index]
    xlim <- CELL_META$xlim
    ylim <- CELL_META$ylim
    xleft <- xlim[1]
    xright <- xlim[2]
    ybottom <- ifelse(strand(meta) == "+", 0, -1)
    ytop <- ifelse(strand(meta) == "+", 1, 0)
    lab <- ifelse(CELL_META$sector.index == "DLP", "CR", CELL_META$sector.index)
    color <- meta$itemRgb
    circos.rect(xleft, ybottom, xright, ytop, col=color)
    if (meta$region %in% c("rRNA", "coding", "D-loop")) {
      textloc <- ifelse(strand(meta) == "+", 0.5, -0.5)
      textcex <- ifelse(meta$name == "MT-ATP8", 0.55, 
                        ifelse(meta$name %in% c("MT-ND3", "MT-ND4L"), 0.7,
                               ifelse(meta$name %in% c("MT-ND6","MT-ATP6"), 0.8,
                                      ifelse(meta$name %in% 
                                             c("MT-RNR2", "MT-ND1", "MT-ND2", 
                                               "MT-CO1", "MT-ND4", "MT-ND5", 
                                               "MT-CYB"), 1.05, 0.95))))

      textfacing <- ifelse(meta$name %in% 
                           c("HVR1","HVR2","HVR3","CR","DLP", 
                             "MT-ATP8","MT-ND6","MT-ND3","MT-ND4L"), 
                           "clockwise", "inside")
      circos.text(mean(xlim), textloc, sub("HVR3", "", lab), 
                  cex=textcex, col="black", facing=textfacing, niceFacing=TRUE)
    }
  }
  dat <- data.frame(name=names(anno), start=start(anno), end=end(anno))
  res <- list(anno=dat, pfun=pfun)
  circos.par("clock.wise"=FALSE, start.degree=90, gap.degree=0)
  circos.genomicInitialize(data=dat, plotType=NULL, major.by=16569)
  circos.track(panel.fun=pfun, ylim=c(-1,1), track.height=0.45, bg.border=NA)
  if (!is.null(variants)) {
    warning("Warning: variant plotting support is very sucky right now...")
  }
  invisible(res)
}
