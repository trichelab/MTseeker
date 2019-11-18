#' plot the (putative) functional impact of mutations to ETP genes as SVG
#'
#' Tim Vickers created a beautiful illustration of the mitochondrial electron
#' transport chain, and that's where coding mitochondrial DNA mutations will 
#' usually hit (we aren't plotting the mitoribosome or tRNAs just yet). So why
#' reinvent the wheel (and possibly make it square)? 
#' 
#' @param variants  an MVRanges or MVRangesList
#' @param defColor  default color (#c9eded is standard) 
#' @param verbose   be verbose, for debugging? (FALSE)
#' 
#' @return          invisibly, the temporary file to which the SVG was written
#'
#' @import     xml2
#' @importFrom utils browseURL
#'
#' @examples
#' library(MTseekerData)
#' data(RONKSvariants) 
#' MTcomplex(RONKSvariants$RO_1)
#' 
#' @export 
MTcomplex <- function(variants, defColor="#c9eded", verbose=FALSE){

  data(mtAnno.rCRS)
  SVGfile <- system.file("extdata", "mtComplexes.svg", package="MTseeker")
  
  mtComplexes <- split(mtAnno, mtAnno$complex)
  complexes <- paste0("complex", names(mtComplexes))
  complexColors <- sapply(mtComplexes, function(x) unique(x$itemRgb))
  names(complexColors) <- complexes 
  ol <- findOverlaps(variants, mtComplexes)
  hits <- sort(unique(names(complexColors)[subjectHits(ol)]))
  misses <- setdiff(complexes, hits)
  
  tmp <- tempfile()
  SVG <- read_xml(SVGfile)
  children <- xml_children(SVG) 
  names(children) <- sapply(children, xml_attr, "id")
  affected_nodes <- grep("complex", names(children), value=TRUE) 
  complexColors[misses] <- defColor

  for (cplx in names(complexColors)) {
    newColor <- complexColors[cplx]
    colorName <- toupper(paste0(cplx, "COLOR"))
    if (verbose) message("Swapping ", colorName, " for ", newColor, "...")
    nodes <- affected_nodes[grep(colorName, 
                                 sapply(children[affected_nodes], 
                                        xml_attr, "style"))]
    for (node in nodes) { 
      if (verbose) message("Fixing ", node, "...")
      xml_attr(children[[node]], "style") <- 
        sub(colorName, newColor, xml_attr(children[[node]], "style"))
      if (verbose) {
        message(node, " style is now ", xml_attr(children[[node]], "style"))
      }
    }
  }

  if (length(hits) > 0) { 
    theGenes <- unique(unlist(mtComplexes[subjectHits(ol)], use.names=FALSE))
    who <- strsplit(as.character(unique(sampleNames(variants))), "\\.")[[1]][1]
    narrative <- paste(length(theGenes), "mitochondrial genic variants in", who)
    xml_text(children[["narrative"]]) <- narrative
  } 

  write_xml(SVG, tmp) 
  browseURL(tmp) 
  invisible(tmp)

}
