#' Locates the variant and determines local bounds within gene
#' Primarily meant to be used within injectMTvarients
#' 
#' @param query      A single variant in MVRanges form 
#' @param coding     Look only at coding regions? (TRUE) 
#'
#' @return           An MVRanges with start/end codons (overlaps return copies)
#'
#' @export
locateMTvariants <- function(query, coding=TRUE) {

  # Checks if the function has been run before
  if ("gene" %in% names(mcols(query)) &
      "region" %in% names(mcols(query)) &
      "localEnd" %in% names(mcols(query)) & 
      "localStart" %in% names(mcols(query)) &
      "startCodon " %in% names(mcols(query)) &
      "endCodon" %in% names(mcols(query))) {
    return(query) # done 
  }
  
  # No input was given
  if (length(query) == 0) return(NULL)
  
  #run serially
  if (is(query, "MVRangesList")) {
    mvr <- MVRangesList(lapply(query, locateVariants))
    return(mvr)
  }
  
  # For now only support rCRS
  stopifnot(genome(query) == "rCRS")
  data("mtAnno.rCRS", package="MTseeker")
  metadata(query)$annotation <- mtAnno
  
  # Localized genic coordinates
  # Ben said he was interested in looking at tRNA regions 
  ### Figure this out later
  #anno <- subset(mtAnno, region %in% c("tRNA", "coding"))

  # decomposeAndCalc throws me an error when I try to use tRNA
  if (coding) anno <- subset(mtAnno, region == "coding")
  else anno <- mtAnno
  ol <- findOverlaps(query, anno, ignore.strand=TRUE)
  
  if (length(ol) == 0) {
    message("No overlapping genes in coding regions for variant: ", names(query))
    return(query[0])
  }
  
  if (length(subjectHits(ol)) > 2) {
    stop("More than 2 overlapping genes in locateVariants")
  }
  
  query$gene <- NA_character_
  query$overlapGene <- NA_character_

  # If there are multiple genes the variant is located within
  # Create multiple copies of the variant to handle each gene that it is located within
  if (length(subjectHits(ol)) > 1) {

    query <- rep(query, length(subjectHits(ol)))
    
    # Assume that there are only 2 overlapping genes
    query[1]$gene <- names(anno)[subjectHits(ol)][1]
    query[2]$gene <- names(anno)[subjectHits(ol)][2]
    query$overlapGene <- paste(names(anno)[subjectHits(ol)], collapse = ",")
  } 
  
  # Otherwise there is no overlap
  else {
    query$gene <- names(anno)[subjectHits(ol)]
  }
  
  # Initialize
  query$region <- NA_character_
  query$localStart <- NA_integer_
  query$localEnd <- NA_integer_
  query$startCodon <- NA_integer_
  query$endCodon <- NA_integer_
  
  for (i in 1:length(query)) {
    
    subHit <- subjectHits(ol)[i]
    
    # Does the variant fall into a coding, tRNA, rRNA, or D-loop region
    query[i]$region <- anno[subHit]$region
    
    # Local start and end (relative to gene of interest)
    # 1-indexed
    query[i]$localStart <- start(query[i]) - start(anno[subHit]) + 1
    query[i]$localEnd <- end(query[i]) - start(anno[subHit]) + 1
    
    ######################### 
    #if (genome(query) == "rCRS" && names(anno[subHit]) == "MT-ND6") {
    #  
    #  subir <- IRanges(query[i]$localStart, query[i]$localEnd)
    #  data(rCRSeq, package="MTseeker")
    #  anno$refDNA <- getSeq(rCRSeq, anno)
    #  refDNA <- anno["MT-ND6"]$refDNA
    #  refSeq <- as.character(unlist(unname(extractAt(anno["MT-ND6"]$refDNA, subir))))
    #  
    #  show(refSeq)
    #  show(ref(query[i]))
    #  
    #  
    #  query[i]$localStart <- query[i]$localStart + 54
    #  query[i]$localEnd <- query[i]$localEnd + 54
    #  
    #  subir <- IRanges(query[i]$localStart, query[i]$localEnd)
    #  refSeq <- as.character(unlist(unname(extractAt(anno["MT-ND6"]$refDNA, subir))))
    #  
    #  show(refSeq)
    #  show(ref(query[i]))
    #  
    #}
    #########################
    
    # Local codon starts and ends
    # Taking a conservative approach
    query[i]$startCodon <- floor(query[i]$localStart / 3)
    query[i]$endCodon <- ceiling(query[i]$localEnd / 3)

    # Everything is 1-indexed (if variant is in the first codon)
    if (query[i]$startCodon == 0) query[i]$startCodon <- 1
    if (query[i]$endCodon == 0) query[i]$endCodon <- 1
  }

  return(query)
}
