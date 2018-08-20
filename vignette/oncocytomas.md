# MTseeker: mitochondrial variant calling, visualization, and annotation in R 

MTseeker works best when given some interesting mitochondria to work with.

First we need to load the oncocytoma BAMs. We don't actually do this here, 
since they are several gigabytes apiece, but notice that all of them have 
been aligned with BWA against the canonical rCRS mitogenome by splicing it 
into hg19.  As opposed to GRCh37, which is what we should have done.

```{r loadBams, eval=FALSE} 
library(MTseeker)

# we use SamBlaster... a lot... in my lab.
BAMfiles <- grep("(split|disc)", value=T, invert=T, list.files(patt=".bam$"))
names(BAMfiles) <- sapply(strsplit(BAMfiles, "\\."), `[`, 1)
BAMs <- data.frame(BAM=BAMfiles, 
                   Sample_Group=ifelse(grepl("NKS",BAMfiles), "normal","tumor"))
rownames(BAMs) <- sub("NKS", "normal", sub("RO", "oncocytoma", rownames(BAMs)))
BAMs$subject <- as.integer(sapply(strsplit(BAMs$BAM, "(_|\\.)"), `[`, 2))

# we merged all the BAMs after-the-fact, so...
BAMs <- subset(BAMs, grepl("merged", BAMs$BAM))
BAMs <- BAMs[order(BAMs$subject), ]

library(parallel) 
options("mc.cores"=detectCores())
MTreads <- getMT(BAMs, filter=FALSE) 
names(MTreads) <- sapply(strsplit(fileName(MTreads), "\\."), `[`, 1)
# saveRDS(MTreads, file="oncocytoma_and_matched_normal_MTreads.rds")

```

We'd like to compute the relative mitochondrial copy number for each:

```{r computeCN, eval=FALSE}
MTreads <- readRDS("MTreads.rds")
mVn <- Summary(MTreads)$mitoVsNuclear
names(mVn) <- names(MTreads) 
CN <- mVn[seq(2,22,2)]/mVn[seq(1,21,2)] 
mtCN <- data.frame(subject=names(CN), CN=CN)

library(ggplot2) 
library(ggthemes)
p <- ggplot(mtCN, aes(x=subject, y=CN, fill=subject)) + 
       geom_col() + theme_tufte() + ylim(0,5) + 
       ylab("Tumor/normal mitochondrial ratio") + 
       ggtitle("Mitochondrial retention in oncocytomas")
ggsave("mtCN.png") 
```

Obviously it's not much good to have a variant caller that can't call variants.

```{r callVariants, eval=FALSE} 
variants <- callMT(MTreads)
saveRDS(variants, file="oncocytoma_and_matched_normal_variants.rds")
```

Show off the results:

```{r plotVariants} 
data(variants, package="MTseeker")
plot(filt(variants))
```

Now let's have some fun:

```{r makeSVG} 
data(variants, package="MTseeker")
SVG <- mtComplex(variants[[2]]) 
library(rsvg) 
rsvg_pdf("RO_1.functionalAnnot.pdf") 
```

And there you have it.
