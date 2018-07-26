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

# add electron transport chain complex annnotations for plotting later
mcols(mtAnno)[, "complex"] <- NA
complexes <- list(
  "I"=c("MT-ND1","MT-ND2","MT-ND3","MT-ND4L","MT-ND4","MT-ND5","MT-ND6"),
  "III"=c("MT-CYB"),
  "IV"=c("MT-CO1","MT-CO2","MT-CO3"),
  "V"=c("MT-ATP8","MT-ATP6")
)
for (i in names(complexes)) mcols(mtAnno)[complexes[[i]], "complex"] <- i

# save the result so we don't have to do it again
save(mtAnno, file="../../data/mtAnno.rCRS.rda")
