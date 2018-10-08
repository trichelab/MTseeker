#' plot a canonical human (or, in principle, any) mitochondrial genome 
#'
#' The default font sizes, orientations, etc. are optimized for a "cold" start;
#' if you want to fiddle with the details, crack open the code and modify it...
#' or alternatively, add sectors/dendrograms inside of this "framed" version.
#'
#' FIXME: add variant type coloration (del=blue, SNV=black, ins=red) 
#' 
#' @param variants  optional MVRanges or MVRangesList to split by strand & plot
#' @param outside   optional MVRanges or MVRangesList to plot outside the circle
#' @param inside    optional MVRanges or MVRangesList to plot inside the circle
#' @param outcol    optional color assignment function or matrix for outside
#' @param incol     optional color assignment function or matrix for inside
#' @param anno      a GRanges (optional, defaults to mtAnno.rCRS if none given)
#' @param how       optional specification for how to plot multiple samples
#' @param ...       other arguments to pass on to called functions
#' 
#' @return          invisibly, a list: `anno` (data.frame) + `pfun` (panel.fun)
#'
#' @import VariantTools
#' @import rtracklayer
#' @import circlize
#' @import viridis
#'
#' @importFrom grDevices col2rgb rgb
#'
#' @examples 
#' library(MTseekerData)
#' data(RONKSvariants) 
#' mtCircos(RONKSvariants)
#' # same as plot(RONKSvariants)
#' title("Renal oncocytomas and normal kidney samples")
#' 
#' @export 
mtCircos <- function(variants=NULL, outside=NULL, inside=NULL, outcol=NULL, 
                     incol=NULL, anno=NULL, how=c("matrix","VAF"), ...) {

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
  circos.par("clock.wise"=FALSE, "start.degree"=90, "gap.degree"=0, 
             "track.margin"=c(0.005, 0.005), "cell.padding"=c(0.005,0,0.005,0), 
             "points.overflow.warning"=FALSE)
  circos.genomicInitialize(data=dat, plotType=NULL, major.by=16569)

  if (!is.null(variants)) {
    message("Splitting variants by strand...")
    stranded <- byStrand(variants)
    message("Replacing `outside` with heavy-strand variants...")
    outside <- stranded$heavy
    message("Replacing `inside` with light-strand variants...")
    inside <- stranded$light
  } 

  # outside track: variant annotations (use granges() if outside is an MVRL)
  if (!is.null(outside)) {
    bed1 <- .makeBed(outside)
    if (is.null(outcol)) outcol <- .newsprint
    circos.genomicHeatmap(bed1, outcol, line_col=.colorCode(bed1$chr), 
                          track.margin=c(0,0), side="outside", border=NA,
                          line_lwd=2) # how to color-code the "matrix" itself?
  } else { 
    circos.track(track.height=0.15, ylim=c(0,1), bg.border=NA)
  }
  
  # main track, gene names and such
  circos.track(panel.fun=pfun, ylim=c(-1,1), track.height=0.5, 
               track.margin=c(0,0), bg.border=NA)

  # inside track: 
  if (!is.null(inside)) {
    bed2 <- .makeBed(inside)
    if (is.null(incol)) incol <- .newsprint
    circos.genomicHeatmap(bed2, incol, line_col=.colorCode(bed2$chr),
                          track.margin=c(0,0), side="inside", border=NA)
  } else { 
    circos.track(track.height=0.15, ylim=c(0,1), bg.border=NA)
  }

  res <- list(anno=dat, pfun=pfun)
  invisible(res)
}




# helper fn
.height <- function(gr) ifelse(gr$region == "tRNA", 0.5, 1)

# helper fn
.halfheight <- function(gr) gr$region == "tRNA"

# helper fn
.textloc <- function(gr) {
  ifelse(.halfheight(gr), .25, .5) * ifelse(strand(gr) == "+", 1, -1)
}

# helper fn
.textbold <- function(gr) ifelse(gr$region %in% c("coding", "rRNA"), 3, 1)

# helper fn
.textcex <- function(gr) { 
  ifelse(gr$name %in% c("HVR1","HVR2","MT-ATP8"), .65,
         ifelse(gr$name %in% c("MT-ND3","MT-ND4L","MT-ND6","MT-CO2","MT-CO3",
                               "MT-RNR1","MT-RNR2","MT-ATP8","MT-ATP6"),.8,.95))
}

# helper fn
.makeBed <- function(x) {
  bed <- switch(class(x),
                "MVRanges"=.mvrToBed(x),
                "MVRangesList"=.mvrlToBed(x),
                "GRanges"=.grToBed(x))
  names(bed)[1] <- "chr"
  return(bed)
}

# helper fn
.mvrToBed <- function(mvr) { 
  bed <- as.data.frame(locateVariants(mvr))[, c("gene", "start", "end")]
  bed$value <- mvr$VAF
  bed <- subset(bed, !is.na(bed[,1]))
  return(bed)
}

# helper fn
.mvrlToBed <- function(mvrl) { 
  gr <- granges(mvrl)
  bed <- as.data.frame(gr)
  bed <- bed[, c(6, 2, 3, 8:(ncol(bed)))]
  return(bed)
}

# helper fn
.grToBed <- function(gr) { 
  bed <- as.data.frame(gr)[, c("gene", "start", "end")]
  bed$value <- 1
  return(bed)
}

# helper fn
.colorCode <- function(x, darken=TRUE, howMuch=1.25) { 
  data("mtAnno.rCRS", package="MTseeker")
  color <- mtAnno[x]$itemRgb
  if (darken) color <- .darken(color, howMuch=howMuch)
  return(color)
}

# helper fn
.darken <- function(hex, howMuch=1.25) {
  rgb(t(col2rgb(hex)/howMuch), maxColorValue=255)
}

# helper fn
.newsprint <- colorRamp2(c(0, 1), c("#FFFFFF", "#000000"))

# helper fn
.bloody <- colorRamp2(c(0, 1), c("#FFFFFF", "#880000"))

# helper fn
.blurple <- colorRamp2(c(0, 1), c("#FFFFFF", "#FF00FF"))

# helper fn
.viridis <- colorRamp2(seq(0, 1, by = 0.1), viridis(11))

# helper fn
.plasma <- colorRamp2(seq(0, 1, by = 0.1), plasma(11))

# helper fn; should probably use viridis instead 
.jet <- colorRamp2(seq(0, 1, 0.125),
                   c("#00007F", "blue", "#007FFF", "cyan",
                     "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

