#' build and install GmapGenome.Hsapiens.[mtGenome], currently for rCRS
#'
#' gmapR needs a reference genome index against which to call variants.
#' We support (only) rCRS as that reference (at least for the time being). 
#' This function creates and installs a reference Gmap index, by default rCRS.
#' In principle, hg19 and mm10 are supported, but in practice, support is poor.
#'
#' @param  mtGenome   mitochondrial reference genome to index (default is rCRS)
#' @param  chrM       what the mitochondrial contig is called (default is chrM)
#' @param  destDir    optional destination for the package ($HOME is default)
#' @param  install    install the package after creation? (default is TRUE) 
#'
#' @return            the path to the created package as a character string
#' 
#' @import gmapR
#'
#' @export
indexMtGenome <- function(mtGenome="rCRS", chrM="chrM", 
                          destDir=NULL, install=TRUE) {

  fa <- system.file(paste0("extdata/", mtGenome, ".fa"), package="MTseeker")
  gmapGenomeRef <- GmapGenome(FastaFile(fa), create=TRUE)
  show(gmapGenomeRef)

  homeDir <- Sys.getenv("HOME") 
  if (is.null(destDir)) destDir <- homeDir
  pkgName <- paste0("GmapGenome.Hsapiens.", mtGenome, ".", chrM)
  pkgPath <- paste0(destDir, "/", pkgName)
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
