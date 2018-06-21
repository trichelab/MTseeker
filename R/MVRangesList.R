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
#' @export
MVRangesList <- function(...) {
  new("MVRangesList", GenomicRangesList(...), elementType = "MVRanges")
}


#' MVRangesList methods (centralized).
#'
#' `counts`               returns fragment counts, if any
#' `counts<-`             adds or updates fragment counts
#' `coverage`             returns estimated coverage for each element
#' `genome`               returns the genome (or, perhaps, genomes) in an MVRL
#' `filt`                 removes variants where PASS != TRUE for each element 
#' `encoding`             returns mutations in coding regions for each element
#' `granges`              returns mildly annotated aggregates of variant sites
#' `tallyVariants`        return a matrix of variant types by annotated region
#' `locateVariants`       locates variants within genes, tRNA, rRNA, or D-loop
#' `summarizeVariants`    attempts mass functional annotation of variant sites
#' 
#' @param x             an MVRangesList (for some methods)
#' @param value         a RangedSummarizedExperiment with matching colnames
#' @param query         an MVRangesList (for predictCoding)
#' @param object        an MVRangesList (for other methods)
#' @param filterLowQual opt. for `granges`/`summarizeVariants`/`tallyVariants`
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
            data(chrominfo.rCRS)
            # pull in annotations
            anno <- suppressMessages(getAnnotations(annotation(x[[1]]))) 
            if (filterLowQual == TRUE) x <- filt(x) 
            message("Aggregating variants...")
            gr <- unlist(as(x, "GRangesList")) 
            seqlevelsStyle(gr) <- "UCSC" # chrM
            gr <- keepSeqlevels(gr, 
                                seqlevels(chrominfo.rCRS), 
                                pruning.mode="coarse")
            mtGenome <- unique(genome(gr))
            gr <- reduce(gr)
            annoGenome <- unique(genome(anno))
            newMtGenome <- unique(genome(gr))
            stopifnot(newMtGenome == annoGenome)
            metadata(gr)$annotation <- anno
            ol <- findOverlaps(gr, anno)
            message("Annotating variants by region...")
            gr$gene <- NA_character_
            gr[queryHits(ol)]$gene <- names(anno)[subjectHits(ol)] 
            gr$region <- NA_character_
            gr[queryHits(ol)]$region <- anno[subjectHits(ol)]$region
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
                res$genomic <- with(res, paste0("g.", Start, Ref, ">", Alt))
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
