#' like a VRanges, but for mitochondria
#' 
#' @import VariantAnnotation
#' @importFrom S4Vectors endoapply
#' @importFrom VariantAnnotation filt 
#' 
#' @exportClass MVRanges
setClass("MVRanges", 
         representation(coverage="numeric"),
         contains="VRanges")


#' wrap a VRanges for mitochondrial use
#'
#' Usually the MVRanges constructor will be called by callMT(). 
#' 
#' @rdname    MVRanges-methods
#'
#' @param   vr        the VRanges
#' @param   coverage  estimated coverage
#'
#' @return            an MVRanges
#'
#' @import BiocGenerics
#'
#' @examples
#' 
#' library(MTseekerData)
#' BAMdir <- system.file("extdata", "BAMs", package="MTseekerData")
#' BAMs <- paste0(BAMdir, "/", list.files(BAMdir, pattern="^pt.*bam$"))
#' (mvr <- filterMTvars(pileupMT(BAMs[1], ref="rCRS")))
#' locateVariants(head(encoding(mvr))) # FIXME: no gene overlaps?!?
#' predictCoding(mvr) # FIXME: none!
#' 
#' # summarizeVariants can take a LONG time to run, and requires Internet
#'
#' @export
MVRanges <- function(vr=NULL, coverage=NA_real_) {
  if (is.null(vr)) vr <- VRanges()
  new("MVRanges", vr, coverage=coverage)
}


#' MVRanges methods (centralized).
#'
#' Many of these methods can be dispatched from an MVRangesList OR an MVRanges.
#' In such cases, the method will usually, but not always, be apply()ed. 
#' 
#' @section Utility methods:
#' 
#' `pos` returns a character vector describing variant positions. 
#' `filt` returns a subset of variant calls where PASS == TRUE (i.e. filtered)
#' `coverage` returns an Rle of coverage across the mitochondrial genome
#' `genomeCoverage` returns the estimated mitochondrial read coverage depth
#'
#' @section Annotation methods:
#' 
#' `type` returns a character vector describing variant type (SNV or indel)
#' `genes` retrieves a GRanges of mitochondrial gene locations for an MVRanges
#' `snpCall` retrieves single nucleotide variant polymorphisms PASSing filters
#' `annotation` gets (perhaps oddly) an MVRanges object annotated against rCRS
#' `getAnnotations` returns the GRanges of gene/region annotations for an MVR
#' `encoding` returns variants residing in coding regions (consequence unknown)
#' `locateVariants` annotates variants w/region, gene, and localStart/localEnd
#' `predictCoding` returns variants consequence predictions as one might expect
#' `summarizeVariants` uses MitImpact to attempt annotation of coding variants.
#' `consensusString` edits rCRS to create a consensus genotype for eg Haplogrep
#'
#' @section Visualization methods:
#'
#' `plot` creates a circular plot of mitochondrial variant calls with annotation
#'
#' @param x             an MVRanges
#' @param object        an MVRanges
#' @param annotations   an MVRanges
#' @param query         an MVRanges
#' @param mode          miscellaneous arguments
#' @param y             another MVRanges
#' @param varAllele     variant alleles
#' @param subject       a GRanges, usually 
#' @param seqSource     a BSgenome, usually 
#' @param ...           miscellaneous args, passed through
#' 
#' @param filterLowQual boolean; drop non-PASSing variants from locateVariants?
#'
#' @return              depends on the method invoked.
#' 
#' @aliases locateVariants getAnnotations predictCoding genes
#' @aliases snpCall annotation summarizeVariants
#'
#' @import              Homo.sapiens
#' 
#' @importFrom          GenomicFeatures   genes
#' @importFrom          Biobase           snpCall
#' @importFrom          graphics          plot 
#'
#' @importMethodsFrom   VariantAnnotation filt 
#' @importMethodsFrom   GenomicFeatures   genes
#' @importMethodsFrom   Biostrings        consensusString type 
#' @importMethodsFrom   IRanges           coverage
#' 
#' @name                MVRanges-methods
NULL


#' @rdname    MVRanges-methods
#' @export
setMethod("genomeCoverage", signature(x="MVRanges"), function(x) x@coverage)


#' @rdname    MVRanges-methods
#' @export
setMethod("coverage", signature(x="MVRanges"), 
          function(x) callNextMethod()[[1]]) 


#' @rdname    MVRanges-methods
#' @export
setMethod("type", signature(x="MVRanges"), 
          function(x) ifelse(nchar(ref(x)) == nchar(alt(x)), "SNV", "indel"))


#' @rdname    MVRanges-methods
#' @export
setMethod("genes", signature(x="MVRanges"), 
          function(x) subset(getAnnotations(x), region == "coding"))


#' @rdname    MVRanges-methods
#' @export
setMethod("snpCall", signature(object="MVRanges"),
          function(object) subset(object, nchar(alt) == nchar(ref)))


#' @rdname    MVRanges-methods
#' @export
setMethod("pos", signature(x="MVRanges"), 
          function(x) {
            .getse <- function(y) cbind(start(y), end(y))
            .condense <- function(y) paste(unique(y), collapse="_")
            apply(.getse(granges(x)), 1, .condense)
          })


#' @rdname    MVRanges-methods
#' @export
setMethod("show", signature(object="MVRanges"),
          function(object) {
            callNextMethod()
            cat(paste0("  genome: ", genome(object)))
            if ("annotation" %in% names(metadata(object))) {
              cat(" (try getAnnotations(object))")
            }
            cat(paste0(", ~",round(genomeCoverage(object)),"x read coverage")) 
            cat("\n")
          })


#' @rdname    MVRanges-methods
#' @export
setMethod("annotation", signature(object="MVRanges"), 
          function(object) {
            
            if (!"annotation" %in% names(metadata(object))) {
              data(mtAnno.rCRS)
              metadata(object)$annotation <- mtAnno
            }
            
            anno <- getAnnotations(object)
            ol <- findOverlaps(object, anno)
            object$gene <- NA_character_
            object[queryHits(ol)]$gene <- names(anno)[subjectHits(ol)] 
            object$region <- NA_character_
            object[queryHits(ol)]$region <- anno[subjectHits(ol)]$region
            return(object)
            
          })


# previously defined in chromvar
setGeneric("getAnnotations",
           function(annotations, ...) standardGeneric("getAnnotations"))


#' @rdname    MVRanges-methods
#' @export
setMethod("getAnnotations", signature(annotations="MVRanges"), 
          function(annotations) {
            if (is.null(metadata(annotations)$annotation)) {
              # FIXME: mouse anno
              data(mtAnno.rCRS)
              return(mtAnno)
            } else { 
              return(metadata(annotations)$annotation)
            }
          })


#' @rdname    MVRanges-methods
#' @export
setMethod("encoding", signature(x="MVRanges"), 
          function(x) {
            
            # limit the search 
            anno <- getAnnotations(x)
            if (!is.null(anno)) {
              x <- subsetByOverlaps(x, subset(anno, region == "coding"))
            } else {
              x <- locateVariants(x) 
              x <- subset(x, region == "coding") 
            }
            chrM <- grep("(MT|chrM|rCRS|RSRS)", seqlevelsInUse(x), value=TRUE)
            return(keepSeqlevels(x, chrM, pruning.mode="coarse"))
            
          })


#' @rdname    MVRanges-methods
#' @importFrom VariantAnnotation filt
#' @export
setMethod("filt", signature(x="MVRanges"), function(x) subset(x, x$PASS==TRUE))


#' @rdname    MVRanges-methods
#' @export
setMethod("genome", signature(x="MVRanges"), 
          function(x) unique(seqinfo(x)@genome))


#' @rdname    MVRanges-methods
#' @export
setMethod("locateVariants", 
          signature(query="MVRanges","missing","missing"),
          function(query, filterLowQual=FALSE, ...) {
           
            if (filterLowQual == TRUE) query <- filt(query)
            if (length(query) == 0) {
              return(query)
            } else if (length(query) > 1) {
              res <- MVRanges()  
              # lapply/getListElement misbehave
              for (i in seq_along(query)) {
                res <- c(res, locateMTvariants(query[i], ...))         
              } 
              return(res)
            } else {
              locateMTvariants(query, ...)
            }
            
          })


#' @rdname    MVRanges-methods
#' @export
setMethod("predictCoding", # mitochondrial annotations kept internally
          signature(query="MVRanges", "missing", "missing", "missing"), 
          function(query, ...) injectMTVariants(filt(query)))


#' @rdname    MVRanges-methods
#' @import    jsonlite
#' @export
setMethod("summarizeVariants", signature(query="MVRanges","missing","missing"),
          function(query, ...) {
            
            # helper function  
            getImpact <- function(pos) {
              url <- paste("http://mitimpact.css-mendel.it", "api", "v2.0",
                           "genomic_position", sub("^(g|m)\\.","",pos), sep="/")
              res <- as.data.frame(read_json(url, simplifyVector=TRUE)$variants)
              if (nrow(res) > 0) {
                res$genomic <- with(res, paste0("m.", Start, Ref, ">", Alt))
                res$protein <- with(res, paste0("p.",AA_ref,AA_position,AA_alt))
                res$change <- with(res, paste(Gene_symbol, protein))
                res[, c("genomic","protein","change","APOGEE_boost_consensus",
                        "MToolBox","Mitomap_Phenotype","Mitomap_Status",
                        "OXPHOS_complex","dbSNP_150_id","Codon_substitution")]
              } else {
                return(NULL)
              }
            }
            
            hits <- lapply(pos(encoding(query)), getImpact)
            hits <- hits[which(sapply(hits, length) > 0)] 
            
            # be precise, if possible
            for (h in names(hits)) {
              j <- hits[[h]]
              if (h %in% j$genomic) hits[[h]] <- j[which(j$genomic == h),]
            }
            
            do.call(rbind, hits)
            
          })


#' @rdname    MVRanges-methods
#' @export
setMethod("plot", signature(x="MVRanges"), 
          function(x, ...) MTcircos(x, ...))


#' @rdname    MVRanges-methods
#' @export
setMethod("consensusString", signature(x="MVRanges"), 
          function(x, ...) {
            supported <- c("rCRS")
            actual <- unique(genome(x))
            stopifnot(unique(genome(x)) %in% supported)
            data(rCRSeq, package="MTseeker")
            mvr <- snpCall(filt(x)) # gross
            alts <- DNAStringSet(replaceAt(rCRSeq[[1]], ranges(mvr), alt(mvr)))
            names(alts) <- actual # genome 
            return(alts)
          })
