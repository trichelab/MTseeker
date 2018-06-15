#' Set gene names for bp ranges. Original code by Stephen Turner at UVA. 
#'
#' FIXME: this is only for hg19/GRCh37 at the moment, need to add hg38/GRCh38
#' 
#' @param mtCalls   a VRanges of mitochondrial calls with VAF and PASS mcols
#' 
#' @return          a coordinate-sorted VRanges with additional annotations
#' 
#' @import VariantTools
#'
#' @export 
addMtGeneLabel <- function(mtCalls) {

  message("Warning: addMtGeneLabel is outdated.")
  message("Use built-in rCRS annotations instead!")

  mtCalls$bp <- start(mtCalls)
  mtCalls$gene <- "Non-Coding" # default 
  mtChr <- as.character(unique(seqnames(mtCalls)))
  mtGenome <- strsplit(unique(na.omit(genome(mtCalls))), "\\.")[[1]][1]
  if (mtGenome == "GRCh37") mtGenome <- "hg19"
  if (mtGenome == "GRCh38") mtGenome <- "hg38"
  data(mtGenes.hg19)
  data(mtGenes.hg38)
  mtGenes <- switch(mtGenome, hg19=mtGenes.hg19, hg38=mtGenes.hg38)
  seqinfo(mtGenes) <- seqinfo(mtCalls)[mtChr]
  olaps <- findOverlaps(mtCalls, mtGenes, type="within")
  mcols(mtCalls)[queryHits(olaps), "gene"] <- 
    mcols(mtGenes)[subjectHits(olaps), "name"]
  attr(mtCalls, "mtGenes") <- mtGenes
  return(sort(mtCalls)) 

}
