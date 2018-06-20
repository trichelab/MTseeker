#' build and install GmapGenome.[organism].[mtGenome], mostly for rCRS
#'
#' gmapR needs a reference genome index against which to call variants.
#' We support (only) rCRS as that reference (at least for the time being). 
#' This function creates and installs a reference Gmap index, by default rCRS.
#' In principle, hg19 and mm10 are supported, but in practice, support is poor.
#' Note: this function creates a "skeleton key" rCRS index for contigs named
#' 'chrM', 'MT', 'rCRS', 'NC_012920.1', and/or 'gi|251831106|ref|NC_012920.1|'.
#'
#' @param  mtGenome   mitochondrial reference genome to index (default is rCRS)
#' @param  fa         FASTA file (default is to find included `mtGenome`.fa)
#' @param  organism   organism whose mitochondrial genome is indexed (Hsapiens) 
#' @param  destDir    optional destination for the package ($HOME is default)
#' @param  install    install the package after creation? (default is TRUE) 
#'
#' @return            the path to the created package as a character string
#' 
#' @import gmapR
#'
#' @export
indexMtGenome <- function(mtGenome="rCRS", fa=NULL, organism="Hsapiens", 
                          destDir=NULL, install=TRUE){

  if (is.null(fa)) { 
    fa <- system.file(paste0("extdata/", mtGenome, ".fa"), package="MTseeker")
  }
  faf <- FastaFile(fa)
  message("Found contigs named ", paste(seqlevels(faf), collapse=", "), "...")
  gmapGenomeRef <- GmapGenome(faf, create=TRUE)
  show(gmapGenomeRef)

  if (is.null(destDir)) destDir <- Sys.getenv("HOME") 
  pkgName <- paste("GmapGenome", organism, mtGenome, sep=".")
  pkgPath <- paste(destDir, pkgName, sep="/")
  message("Building ", pkgName, " in ", pkgPath)
  makeGmapGenomePackage(gmapGenomeRef, 
                        version="1.0", 
                        author="Me", 
                        maintainer="My Self <myself@me.com>", 
                        destDir=destDir, 
                        pkgName=pkgName, 
                        license="Artistic-2.0")
  if (install == TRUE) install.packages(pkgPath, repos=NULL)
  return(pkgPath)

}
