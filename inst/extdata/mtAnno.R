library(MTseeker)
library(rtracklayer)
mtColors <- read.csv(system.file("extdata", "mtColors.csv", package="MTseeker"),
                     row.names=1, stringsAsFactors=FALSE) 

# turn it into a BED9
mtAnno <- makeGRangesFromDataFrame(subset(mtColors, !is.na(chrom)), keep=TRUE)
names(mcols(mtAnno)) <- sub("color", "itemRgb", names(mcols(mtAnno))) 
export(mtAnno, "mtAnno.rCRS.bed")

# turn into annotations!
genome(mtAnno) <- "rCRS"
seqlengths(mtAnno) <- 16569
isCircular(mtAnno) <- TRUE
save(mtAnno, file="../../data/mtAnno.rCRS.rda")
