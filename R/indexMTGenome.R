#' build and install GmapGenome.[organism].[mtGenome], currently Hsapiens.rCRS
#'
#' gmapR needs a reference genome index in order for it to call any variants.
#' We support rCRS (and only rCRS) as that reference, at least for the moment.
#' This function creates & installs a reference (rCRS, the default) Gmap index.
#' In principle, hg19 and mm10 could be supported; in practice, support is poor.
#' (Also, the Yoruban chrM in hg19 is a terrible reference for variant calling.)
#' We would be grateful for a patch to add mm10/GRCm38 support; eventually, we 
#' plan to add it in ourselves (as one might have guessed from Mouse Mitocarta).
#' 
#' Note: this function creates a "skeleton key" rCRS index for contigs named
#' 'chrM', 'MT', 'rCRS', 'NC_012920.1', and/or 'gi|251831106|ref|NC_012920.1|'.
#' The point of this kludge is to allow gmapR to call variants against various
#' styles of contig names, whether NCBI, UCSC, Genbank, or colloquial rCRS.
#'
#' @param  mtGenome   mitochondrial reference genome to index (default is rCRS)
#' @param  fa         FASTA file (default is to find included `mtGenome`.fa)
#' @param  organism   organism whose mitochondrial genome is indexed (Hsapiens) 
#' @param  destDir    optional destination for the package ($HOME is default)
#' @param  install    install the package after creation? (default is TRUE) 
#' @param  unlink     if an index package already exists, remove it? (FALSE)
#'
#' @return            the path to the created package as a character string
#'
#' @examples 
#' 
#' if (.Platform$OS.type != "windows") {
#'   mtGenome <- "rCRS"
#'   fa <- system.file(paste0("extdata/", mtGenome, ".fa"), package="MTseeker")
#'   indexMTGenome(mtGenome=mtGenome, fa=fa, destDir=tempdir())
#' }
#'
#' @import gmapR
#' @importFrom utils install.packages installed.packages remove.packages
#'
#' @export
indexMTGenome <- function(mtGenome="rCRS", fa=NULL, organism="Hsapiens", 
                          destDir=NULL, install=TRUE, unlink=FALSE) {

  # default to rCRS
  if (is.null(fa)) { 
    fa <- system.file(paste0("extdata/", mtGenome, ".fa"), package="MTseeker")
  }
  faf <- FastaFile(fa)
  fai <- scanFaIndex(fa)
  message("Found contigs named ", paste(seqlevels(fai), collapse=", "), "...")

  # stick it in $HOME if no better option is provided
  if (is.null(destDir)) destDir <- Sys.getenv("HOME") 
  pkgName <- paste("GmapGenome", organism, mtGenome, sep=".")
  pkgPath <- paste(destDir, pkgName, sep="/")
  cachePath <- paste(destDir, ".local/share/gmap")

  # replace? 
  if (unlink) { 
    if (dir.exists(pkgPath)) unlink(pkgPath, recursive=TRUE) 
    if (dir.exists(cachePath)) unlink(cachePath, recursive=TRUE) 
    if (pkgName %in% rownames(installed.packages()) & install) {
      message("Removing installed ", pkgName, "...") 
      unloadNamespace(pkgName)
      remove.packages(pkgName)
    }
  }

  # index the provided FASTA file 
  gmapGenomeRef <- GmapGenome(faf, create=TRUE)
  show(gmapGenomeRef)

  # the author/maintainer fields are a kludge
  message("Building ", pkgName, " in ", pkgPath)
  makeGmapGenomePackage(gmapGenomeRef, 
                        version="1.0", 
                        author="Me", 
                        maintainer="My Self <myself@me.com>", 
                        destDir=destDir, 
                        pkgName=pkgName, 
                        license="Artistic-2.0")

  # this raises a NOTE: in BiocCheck; I think it's worth it. 
  if (install == TRUE) install.packages(pkgPath, repos=NULL)
  return(pkgPath)

}
