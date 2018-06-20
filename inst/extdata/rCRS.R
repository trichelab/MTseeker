library(gmapR)
library(MTseeker)
library(rtracklayer)

fa <- system.file("extdata/rCRS.fa", package="MTseeker")
gmapGenomeRCRS <- GmapGenome(FastaFile(fa), create=TRUE)
show(gmapGenomeRCRS)

makeGmapGenomePackage(gmapGenomeRCRS, 
                      version="1.0", 
                      author="Me", 
                      maintainer="My Self <myself@me.com>", 
                      destDir=".", pkgName="GmapGenome.Hsapiens.rCRS",
                      license="Artistic-2.0")
