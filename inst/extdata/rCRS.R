library(gmapR)
library(ATACseeker)
library(rtracklayer)

fa <- system.file("extdata/mitomes/rCRS.fasta", package="ATACseeker")
gmapGenomeRCRS <- GmapGenome(FastaFile(fa), create=TRUE)
show(gmapGenomeRCRS)

makeGmapGenomePackage(gmapGenomeRCRS, 
                      version="0.99", 
                      author="Me", 
                      maintainer="My Self <myself@me.com>", 
                      destDir=".", pkgName="GmapGenome.Hsapiens.rCRS",
                      license="Artistic-2.0")
