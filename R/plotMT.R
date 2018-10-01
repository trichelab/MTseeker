#' plot mitochondrial variant calls in the context of heteroplasmy
#'
#' use a plot (inspired by the original, which was written by Stephen Turner)
#' to visualize mitochondrial heteroplasmy from mitochondrial variant calls
#'
#' by default, calls with a VAF of 0 or 1 are excluded, but this can be changed
#'
#' @param mtCalls     a VRanges object with variant calls and associated VAFs
#' @param filterVAF   filter low-quality and 0 or 1 VAF calls? (TRUE)
#' @param rot         plot counterclockwise (-1, default) or clockwise (1)?
#' @param title       a title for the plot (NULL)
#' 
#' @return a ggplot object
#' 
#' @import ggplot2
#' @import ggthemes
#' @import VariantTools
#' 
#' @export
plotMT <- function(mtCalls, filterVAF=TRUE, rot=-1, title=NULL) {

  if (!is(mtCalls, "VRanges")) stop("plotMT() expects an MVRanges object.")
  if (filterVAF == TRUE) {
    mtCalls <- subset(mtCalls, 
                      mtCalls$VAF!=1 & mtCalls$VAF!=0 & mtCalls$PASS==1)
  }
  mtCalls <- addMtGeneLabel(mtCalls)
  mtGenes <- attr(mtCalls, "mtGenes") # for boundaries
  mtChr <- seqlevels(mtGenes)[1]
  stopifnot(is(mtGenes, "GRanges"))
  stopifnot(rot %in% c(-1, 1))

  # Set colors for each gene
  colours <- c("Control-Region"="lightblue4",
               "tRNA"="magenta4", 
               "rRNA"="mediumaquamarine",
               "Non-Coding"="sienna4", 
               "MT-ND1"="magenta", 
               "MT-ND2"="mediumblue", 
               "MT-CO1"="olivedrab", 
               "MT-CO2"="orange2", 
               "MT-ATP8"="orchid4", 
               "MT-ATP6"="red3", 
               "MT-CO3"="royalblue2", 
               "MT-ND3"="palegreen4",
               "MT-ND4L"="grey0", 
               "MT-ND4"="pink4", 
               "MT-ND5"="yellow4",
               "MT-ND6"="steelblue4", 
               "MT-CYB"="tan",
               "red")

  # Create gene boundaries and lines
  mtSeqLen <- seqlengths(mtGenes)[1]
  mtRegions <- subset(mtGenes, !mtGenes$name %in% c("tRNA", "rRNA"))
  visibleboundaries <- start(mtRegions) + ((end(mtRegions)-start(mtRegions))/2)
  bdries <- data.frame(bp=visibleboundaries, VAF=-0.25)
  bdries$gene <- mtRegions$name
  bdries$angle <- with(bdries, (90 + (360*bp/mtSeqLen)) %% 360)

  gap <- gaps(mtGenes)
  mtNonCoding <- subset(gap, decode(strand(gap)) == "*")
  mtNonCoding$name <- "Non-Coding"
  mtAll <- sort(c(mtGenes, mtNonCoding))

  l <- data.frame(bp=seq(0, mtSeqLen), VAF=0)
  l$gene <- c("Control-Region", 
              unlist(sapply(mtAll, function(x) rep(x$name, width(x)))))

  colourBreaks <- names(colours)[-length(colours)]
  colourLabels <- colourBreaks

  whatKindOf <- ifelse(filterVAF, "Heteroplasmic", "Mitochondrial")
  plt <- ggplot(as.data.frame(mcols(mtCalls)), 
                aes(x=bp, xend=bp, y=VAF, yend=0, color=gene)) +
         theme_tufte() + 
         geom_segment(size=0.25, show.legend=TRUE) + 
         scale_y_continuous(labels=scales::percent, 
                            breaks=c(0, 0.25, 0.5, 0.75, 1), 
                            limits=c(-3, 1)) + 
         coord_polar(direction=rot) +  # can change this
         scale_color_manual(values=colours,
                            breaks=colourBreaks,
                            labels=colourLabels) +
         geom_text(aes(label=gene, angle=angle), check_overlap=T, 
                   data=bdries, size=3, hjust=1) +
         geom_point(aes(color=gene), data=l) + 
         ggtitle(paste(whatKindOf, "variants: position and VAF")) + 
         guides(colour=guide_legend(title="Gene")) + 
         xlab("") +
         ylab("")  
  
  return(plt)
}
