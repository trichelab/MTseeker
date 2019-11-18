#' @rdname    MTcircos
#'
#' @export
genMTcircos <- function(mvr) {
  
  circos.clear()
  
  anno <- initMTcircos(mvr)
  genesMTcircos(mvr, anno, legends=TRUE)
  
}

#' @rdname    MTcircos
#'
#' @export
initMTcircos <- function(x) {
  
  # Human
  if (unique(genome(x)) == "rCRS") {
    data(mtAnno.rCRS)
    anno <- mtAnno #.rCRS
    refWidth <- 16569
  }
  
  # Mouse
  else if (unique(genome(x)) == "NC_005089") {
    # FIXME: this will always fail until it is moved into MTseeker/data 
    anno <- readRDS("~/Documents/pileupTesting/NC_005089genome/MTmouseAnno.rds")
    refWidth <- 16299
  }
  
  dat <- data.frame(name=names(anno), start=start(anno), end=end(anno))
  
  circos.par("clock.wise"=FALSE, "start.degree"=90, "gap.degree"=0, 
             "track.margin"=c(0.005, 0.005), "cell.padding"=c(0.005,0,0.005,0), 
             "points.overflow.warning"=FALSE)
  
  # Important to set the factor levels otherwise plot is out of order
  circos.genomicInitialize(factor=factor(dat$name, levels=dat$name), 
                           data=dat, plotType=NULL)
  
  return(anno)
}


#' @rdname    MTcircos
#' 
#' @param x       something with a genome() 
#' @param anno    a GRanges with annotations
#' @param legends plot legends? (FALSE) 
#'
#' @return        annotations
#' 
#' @export
genesMTcircos <- function(x, anno, legends=F) {
  
  dat <- data.frame(name=names(anno), 
                    start=start(anno), 
                    end=end(anno), 
                    stringsAsFactors = F)
  
  row.names(dat) <- names(anno)

  pfun <- function(x, y) {
    
    xlim <- CELL_META$xlim
    ylim <- CELL_META$ylim
    
    gr <- anno[CELL_META$sector.index]
    
    ytop <- .height(gr) * ifelse(strand(gr) == "+", 1, 0)
    ybot <- .height(gr) * ifelse(strand(gr) == "-", -1, 0)
    lab <- ifelse(CELL_META$sector.index == "DLP", "CR", CELL_META$sector.index)
    
    circos.rect(xlim[1], ybot, xlim[2], ytop, col=gr$itemRgb)
    
    # Do not assign names to tRNA and HVR3 because these regions are too small
    if (gr$region %in% c("rRNA", "coding", "D-loop") & gr$name != "HVR3") {
      circos.text(mean(xlim), .textloc(gr), lab, col="black", cex=.textcex(gr),
                  font=.textbold(gr), facing="clockwise", niceFacing=TRUE)
    }

  }

  circos.track(panel.fun=pfun, ylim=c(-1,1), track.height=0.5, 
               track.margin=c(0,0), bg.border=NA)
  
  if (legends && genome(x) == "rCRS") {
    
    # Get the colors to correspond which region they belong to
    # Only support rCRS 
    colDF <- as.data.frame(unique(anno$itemRgb), stringsAsFactors=F)
    names(colDF) <- "col"
    
    colDF$label <- NA_character_
    colDF$label[1] <- "Control Region"
    colDF$label[2] <- "D Loop"
    colDF$label[3] <- "tRNA"
    colDF$label[4] <- "rRNA"
    colDF$label[5] <- "Complex I"
    colDF$label[6] <- "Complex IV"
    colDF$label[7] <- "Complex V"
    colDF$label[8] <- "Complex III"
    
    legend("topleft", ncol=2, title="MT Regions",
           legend=colDF$label, col=colDF$col, pch=15, cex=0.6)
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
