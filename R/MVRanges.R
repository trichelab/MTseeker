#' like a VRanges, but for mitochondria
#' 
#' @import VariantAnnotation
#' 
#' @exportClass MVRanges
setClass("MVRanges", 
         representation(coverage="numeric"),
         contains="VRanges")


#' wrap a VRanges for mitochondrial use
#'
#' @param   vr        the VRanges
#' @param   coverage  estimated coverage
#'
#' @return            an MVRanges
#' 
#' @export
MVRanges <- function(vr, coverage) new("MVRanges", vr, coverage=coverage)


#' MVRanges methods (centralized).
#'
#' `pos` returns a character vector describing variant positions 
#' `type` returns a character vector describing variant type (SNV or indel)
#' `filt` returns a subset of variant calls where PASS == TRUE (i.e. filtered)
#' `coverage` returns the estimated average mitochondrial read coverage depth
#' `annotation` returns (perhaps oddly) an annotated, lifted MVRanges object
#' `getAnnotations` returns the GRanges of gene/region annotations for an MVR
#' `encoding` returns variants residing in coding regions (consequence unknown)
#' `predictCoding` returns variants consequence predictions as one might expect
#' `tallyVariants` returns a named vector of variant types by annotated region.
#' `summarizeVariants` uses MitImpact to attempt annotation of coding variants.
#'
#' @param x             an MVRanges
#' @param object        an MVRanges
#' @param annotations   an MVRanges
#' @param query         an MVRanges
#'
#' @name                MVRanges-methods
NULL


#' @rdname    MVRanges-methods
#' @export
setMethod("coverage", signature(x="MVRanges"), function(x) x@coverage)


#' @rdname    MVRanges-methods
#' @export
setMethod("type", signature(x="MVRanges"), 
          function(x) ifelse(width(x) == 1, "SNV", "indel"))


#' @rdname    MVRanges-methods
#' @export
setMethod("pos", signature(x="MVRanges"), 
          function(x) {
            mtChr <- grep("(chrM|MT|rCRS|NC_012920.1)",seqlevels(x),value=TRUE)
            loci <- gsub(paste0(mtChr,":"), "", as.character(granges(x)))
            return(loci)
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
            cat(paste0(", ~", round(coverage(object)), "x read coverage")) 
            cat("\n")
          })


#' @rdname    MVRanges-methods
#' @export
setMethod("annotation", signature(object="MVRanges"), 
          function(object) {

            if (!"annotation" %in% names(metadata(object))) {
              data(anno_rCRS)
              metadata(object)$annotation <- anno_rCRS
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
            anno <- metadata(annotations)$annotation
            if (is.null(anno)) {
              message("Unannotated! Try getAnnotations(annotation(object))).")
            }
            return(anno)
          })


#' @rdname    MVRanges-methods
#' @export
setMethod("encoding", signature(x="MVRanges"), 
          function(x) {

            # limit the search 
            x <- annotation(x) # ensure it's rCRS
            x <- subset(x, x$PASS == TRUE & region == "coding") 

            # fix issues
            data(rCRSeq)
            x <- keepSeqlevels(x, names(rCRSeq), pruning.mode="coarse")
            comp <- data.frame(ref=as.character(getSeq(rCRSeq, x)), alt=alt(x))
            keep <- with(comp, which(ref != alt))
            
            # return subset
            return(x[keep])

          })


#' @rdname    MVRanges-methods
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
          function(query, filterLowQual=TRUE, ...) {

            if (filterLowQual == TRUE) query <- filt(query)
            if ("gene" %in% names(mcols(query)) &
                "region" %in% names(mcols(query))) return(query) # done 

            # otherwise, annotate:
            isCircular(query)["chrM"] <- TRUE # grrr
            whichGenes <- paste0("mtGenes.", genome(query))
            avail <- data(package="MTseeker")$results[,"Item"]
            availableMtGenes <- grep("^mtGenes\\.", value=TRUE, avail)
            if (!whichGenes %in% availableMtGenes) {
              stop("No quick gene location database for ", whichGenes)
            } else { 
              data(list=whichGenes, package="MTseeker")
              anno <- get(whichGenes)
            }

            ol <- findOverlaps(query, anno, ignore.strand=TRUE)
            query$gene <- NA_character_
            query[queryHits(ol)]$gene <- names(anno)[subjectHits(ol)] 
            query$region <- NA_character_
            query[queryHits(ol)]$region <- anno[subjectHits(ol)]$region
            return(query)

          })


#' @rdname    MVRanges-methods
#' @export
setMethod("tallyVariants", signature(x="MVRanges"), 
          function(x, filterLowQual=TRUE, ...) {

            located <- locateVariants(x, filterLowQual=filterLowQual)
            table(located$region)

          })


#' @rdname    MVRanges-methods
#' @export
setMethod("predictCoding", # mitochondrial annotations kept internally
          signature(query="MVRanges", "missing", "missing", "missing"), 
          function(query, ...) {

            # setup:
            data(rCRSeq)
            query <- encoding(query)
            MT_CODE <- getGeneticCode("SGC1")
            mtGenes <- subset(metadata(query)$annotation, region == "coding")
            mcols(mtGenes) <- DataFrame(DNA=getSeq(rCRSeq, mtGenes))
            mtGenes$AA <- translate(mtGenes$DNA, MT_CODE)
            ol <- findOverlaps(query, mtGenes)

            # execution:
            result <- granges(query)
            result$varAllele <- alt(query)
            stop("predictCoding(MVRanges) is not finished...") 

          })


#' @rdname    MVRanges-methods
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
                        "MtoolBox","Mitomap_Phenotype","Mitomap_Status",
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
