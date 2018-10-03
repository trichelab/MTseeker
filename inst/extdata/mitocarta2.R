library(rtracklayer)
library(Homo.sapiens)
library(Mus.musculus)

# hg19 
mitocarta2.hg19 <- import("ftp://ftp.broadinstitute.org/distribution/metabolic/papers/Pagliarini/MitoCarta2.0/Human.MitoCarta2.0.bed")
seqinfo(mitocarta2.hg19) <- seqinfo(Homo.sapiens)[seqlevels(mitocarta2.hg19)]
names(mitocarta2.hg19) <- mapIds(Homo.sapiens, mitocarta2.hg19$name, 
                                 "SYMBOL", "ACCNUM") 
save(mitocarta2.hg19, file="../../data/mitocarta2.hg19.rda")

# mm10 
mitocarta2.mm10 <- import("ftp://ftp.broadinstitute.org/distribution/metabolic/papers/Pagliarini/MitoCarta2.0/Mouse.MitoCarta2.0.bed")
seqinfo(mitocarta2.mm10) <- seqinfo(Mus.musculus)[seqlevels(mitocarta2.mm10)]
names(mitocarta2.mm10) <- mapIds(Mus.musculus, mitocarta2.mm10$name, 
                                 "SYMBOL", "ACCNUM") 
save(mitocarta2.mm10, file="../../data/mitocarta2.mm10.rda")
