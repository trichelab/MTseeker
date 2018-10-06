#' like a VRangesList, but for mitochondria
#' 
#' @import VariantAnnotation
#' @import VariantTools
#' @import Biostrings
#' @import S4Vectors
#' 
#' @exportClass MVRangesList
setClass("MVRangesList", contains="SimpleVRangesList")


#' wrap a VRangesList for mitochondrial use
#'
#' @param ...     the MVRanges elements forming the MVRangesList
#'
#' @return        the MVRangesList
#' 
#' @examples
#' library(MTseekerData)
#' data(RONKSvariants)
#' show(RONKSvariants)
#'
#' @export
MVRangesList <- function(...) {
  new("MVRangesList", GenomicRangesList(...), elementType = "MVRanges")
}


#' MVRangesList methods (centralized).
#'
#' @section Utility methods:
#' 
#' `counts`               returns fragment counts, if any
#' `counts<-`             adds or updates fragment counts
#' `coverage`             returns estimated mitochondrial read coverage depth
#' `filt`                 removes variants where PASS != TRUE for each element 
#'
#' @section Annotation methods:
#'
#' `genes`                returns an annotated GRanges of mitochondrial genes 
#' `getAnnotations`       returns a GRanges of annotated mitochondrial features
#' `genome`               returns the genome (or, perhaps, genomes) in an MVRL
#' `encoding`             returns mutations in coding regions for each element
#' `granges`              returns mildly annotated aggregates of variant sites
#' `snpCall`              retrieves single nucleotide variant polymorphisms 
#' `tallyVariants`        return a matrix of variant types by annotated region
#' `locateVariants`       locates variants within genes, tRNA, rRNA, or D-loop
#' `summarizeVariants`    attempts mass functional annotation of variant sites
#' 
#' @section Visualization methods:
#'
#' `plot`                 creates circular plot of mitochondrial variant calls
#'
#' @param x             an MVRangesList (for some methods)
#' @param value         a RangedSummarizedExperiment with matching colnames
#' @param query         an MVRangesList (for predictCoding)
#' @param object        an MVRangesList (for other methods)
#' @param annotations   an MVRangesList (for getAnnotations)
#' @param filterLowQual opt. for `granges`/`summarizeVariants`/`tallyVariants`
#'
#' @return              the return value depends on the method invoked.
#' 
#' @name  MVRangesList-methods
NULL


#' @rdname    MVRangesList-methods
#' @export
setMethod("coverage", signature(x="MVRangesList"), 
          function(x) sapply(x, coverage))


#' @rdname    MVRangesList-methods
#' @export
setReplaceMethod("counts", 
                 signature(object="MVRangesList", 
                           value="RangedSummarizedExperiment"),
                 function(object, value) {
                   if (!identical(names(object), colnames(value))) {
                     stop("Error: colnames(value) doesn't match names(object)!")
                   } else if (!"counts" %in% names(assays(value))) {
                     stop("Error: value must have an assay named `counts`!")
                   } else {
                     columns <- names(object)
                     metadata(object)$counts <- filterPeaks(value[, columns])
                     return(object)
                   }
                 })


#' @rdname    MVRangesList-methods
#' @export
setMethod("counts", signature(object="MVRangesList"), 
          # it turns out that filtering may be needed on egress:
          function(object) filterPeaks(metadata(object)$counts))


#' @rdname    MVRangesList-methods
#' @export
setMethod("genes", signature(x="MVRangesList"), 
          function(x) subset(getAnnotations(x), region == "coding"))


#' @rdname    MVRangesList-methods
#' @export
setMethod("snpCall", signature(object="MVRangesList"),
          function(object) endoapply(object, snpCall))


#' @rdname    MVRangesList-methods
#' @export
setMethod("getAnnotations", signature(annotations="MVRangesList"), 
          function(annotations) getAnnotations(annotations[[1]]))


#' @rdname    MVRangesList-methods
#' @export
setMethod("encoding", signature(x="MVRangesList"), 
          function(x) MVRangesList(lapply(x, encoding)))


#' @rdname    MVRangesList-methods
#' @export
setMethod("predictCoding", # mitochondrial annotations kept internally
          signature(query="MVRangesList", "missing", "missing", "missing"), 
          function(query, ...) sapply(query, predictCoding))


#' @rdname    MVRangesList-methods
#' @export
setMethod("show", signature(object="MVRangesList"),
          function(object) {
            callNextMethod()
            coverages <- paste0(round(unname(sapply(object, coverage))), "x")
            cat(S4Vectors:::labeledLine("coverage", coverages))
            if ("counts" %in% names(metadata(object))) {
              peaks <- nrow(metadata(object)$counts)
              cat(ifelse("bias" %in% names(rowData(counts(object))),
                  "Bias-corrected ", "Raw "))
              cat("fragment counts at", peaks, "peaks are available from",
                  "counts(object).\n")
            }
          })


# helper
setAs(from="MVRangesList", to="GRangesList",
      function(from) GRangesList(lapply(from, granges)))


#' @rdname    MVRangesList-methods
#' @export
setMethod("filt", signature(x="MVRangesList"),
          function(x) {
            message("Filtering out low-quality calls...")
            MVRangesList(sapply(x, subset, PASS))
          })


#' @rdname    MVRangesList-methods
#' @export
setMethod("granges", signature(x="MVRangesList"),
          function(x, filterLowQual=TRUE) {

            # if cached...
            if (filterLowQual == TRUE & 
                "granges.filtered" %in% names(metadata(x))) {
              return(metadata(x)$granges.filtered)
            } else if (filterLowQual == FALSE & 
                       "granges.unfiltered" %in% names(metadata(x))) {
              return(metadata(x)$granges.unfiltered)
            }

            # if not...
            if (filterLowQual == TRUE) x <- filt(x) 
            anno <- suppressMessages(getAnnotations(x))
            message("Aggregating variants...")
            gr <- unlist(as(x, "GRangesList")) 
            gr <- keepSeqlevels(gr, "chrM", pruning.mode="coarse")
            gr <- reduce(gr)
            ol <- findOverlaps(gr, anno)
            metadata(gr)$annotation <- anno
            metadata(gr)$sampleNames <- names(x)
            message("Annotating variants by region...")
            gr$gene <- NA_character_
            gr[queryHits(ol)]$gene <- names(anno)[subjectHits(ol)] 
            gr$region <- NA_character_
            gr[queryHits(ol)]$region <- anno[subjectHits(ol)]$region
            message("Annotating variants by sample...") 
            hitMat <- matrix(0, ncol=length(x), nrow=length(gr),
                             dimnames=list(NULL, names(x)))
            varHits <- findOverlaps(as(x, "GRangesList"), gr)
            bySample <- split(subjectHits(varHits), queryHits(varHits))
            for (s in seq_along(bySample)) hitMat[bySample[[s]], s] <- 1
            mcols(gr)$overlaps <- hitMat 
            return(gr)
          })


#' @rdname    MVRangesList-methods
#' @export
setMethod("summarizeVariants", 
          signature(query="MVRangesList","missing","missing"),
          function(query, filterLowQual=TRUE, ...) {
            
            # code duplication! refactor
            getRangedImpact <- function(pos) {
              url <- paste("http://mitimpact.css-mendel.it", "api", "v2.0",
                           "genomic_position", pos, sep="/")
              res <- as.data.frame(read_json(url, simplifyVector=TRUE)$variants)
              if (nrow(res) > 0) {
                res$genomic <- with(res, paste0("m.", Start, Ref, ">", Alt))
                res$protein <- with(res, paste0("p.",AA_ref,AA_position,AA_alt))
                res$Consequence <- with(res, 
                                        paste(Gene_symbol, 
                                              paste0(AA_ref,
                                                     AA_position,
                                                     AA_alt)))
                res[, c("genomic","protein","Start",
                        "Ref","Alt","Codon_substitution","dbSNP_150_id",
                        "Mitomap_Phenotype","Mitomap_Status",
                        "Gene_symbol","OXPHOS_complex",
                        "Consequence","APOGEE_boost_consensus","MtoolBox")]
              } else {
                return(NULL)
              }
            }

            gr <- granges(query, filterLowQual=filterLowQual, ...)
            names(gr) <- as.character(gr)
            message("Retrieving functional annotations for variants...")
            hits <- lapply(as.character(ranges(gr)), getRangedImpact)
            rsv <- do.call(rbind, hits[which(sapply(hits, length) > 0)])
            names(rsv) <- sub("Start", "start", names(rsv)) # grrr
            rsv$chrom <- "chrM"
            rsv$end <- rsv$start # FIXME
            res <- makeGRangesFromDataFrame(rsv, keep=TRUE)
            seqinfo(res) <- seqinfo(gr)
            return(res)

          })


#' @rdname    MVRangesList-methods
#' @export
setMethod("genome", signature(x="MVRangesList"),
          function(x) unique(sapply(lapply(x, genome), unique)))


#' @rdname    MVRangesList-methods
#' @export
setMethod("locateVariants", 
          signature(query="MVRangesList","missing","missing"),
          function(query, filterLowQual=TRUE, ...) {

            stop("Don't use this method for now. It has bugs!")
            if (filterLowQual == TRUE) query <- filt(query)
            MVRangesList(lapply(query, locateVariants))

          })


#' @rdname    MVRangesList-methods
#' @export
setMethod("tallyVariants", signature(x="MVRangesList"),
          function(x, filterLowQual=TRUE, ...) {
            
            stop("Don't use this method for now. It has bugs!")
            t(sapply(x, tallyVariants))

          })


#' @rdname    MVRangesList-methods
#' @export
setMethod("plot", signature(x="MVRangesList"),
          function(x, ...) mtCircos(x, ...))
