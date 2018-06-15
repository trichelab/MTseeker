#' liftOver and simplify an hg19 mitochondrial genome (variants) to rCRS 
#'
#' Note: this presupposes that the rCRS sequence (GenBank accession NC_012920.1)
#' is available and gene features provided by GenBank as a GFF3 are available,
#' or can be located. Otherwise this is all rather a waste of time...
#'
#' @param mvr   an MVRanges of variant calls
#' @param quiet suppress warnings and messages? (FALSE)
#'
#' @return      an MVRanges of variant calls, repositioned against rCRS 
#' 
#' @import      rtracklayer
#' @import      Biostrings
#'
#' @export
rCRS <- function(mvr, quiet=FALSE) { 

  mtGenome <- unique(genome(mvr))
  if (mtGenome == "rCRS") return(mvr)

  data(chrominfo.rCRS)
  seqlevelsStyle(mvr) <- "UCSC" # chrM
  if (mtGenome %in% c("GRCh38","hg38")) {
    seqinfo(mvr) <- chrominfo.rCRS # identical to GRCh38/hg38 modulo name
  } else if (mtGenome == "hg19") { 
    data(hg19TorCRS) 
    mvr <- sort(unlist(liftOver(mvr, hg19TorCRS)))
    seqinfo(mvr) <- chrominfo.rCRS # as with TxDB 
    if (quiet == FALSE) { 
      message("Note: you will have better results if you realign against rCRS.")
    }
  } else { 
    message("Your variants appear to have been called against ", mtGenome, ".")
    stop("This function currently supports (only) hg19 and hg38/GRCh38.")
  }

  names(mvr) <- mtHGVS(mvr)
  mvr <- annotation(mvr)
  return(mvr)

}
