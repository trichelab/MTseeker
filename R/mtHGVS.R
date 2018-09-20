#' convert mitochondrial variant calls to HGVS format for naming
#'
#' @param x         an MVRanges (or, in a pinch, a GRanges)
#' @param verbose   be yappy? (default is FALSE) 
#' 
#' @return          proper HGVS names for the MVRanges (or GRanges)
#'
#' @export
mtHGVS <- function(x, verbose=FALSE) { 

  # NCBI recommendation is ASSEMBLY:m.HGVS
  pre <- paste(unique(genome(x)), "m.", sep=":") 

  # for plotting purposes
  if (class(x) == "GRanges") {
    return(paste0(pre, sapply(apply(cbind(start(x), end(x)), 1, unique), 
                              paste, collapse="-")))
  }

  # baseline notation, assuming mostly SNVs 
  nms <- paste0(gsub("\\-", "_", pos(x)), ref(x), ">", alt(x))

  # helper fn
  psub <- function(x, y) { 
    stopifnot(length(x) == length(y))
    for (i in seq_along(x)) y[i] <- sub(paste0("^", x[i]), "", y[i])
    return(y)
  }

  # helper fn
  plural <- function(n, x, suf="s") paste(n, ifelse(n == 1, x, paste0(x, suf)))
  
  # single nucleotide deletions 
  sbdel <- which((nchar(ref(x)) - nchar(alt(x))) == 1)
  sbdels <- length(sbdel) 
  if (sbdels > 0) {
    nms[sbdel] <- paste0(start(x)[sbdel] + 1, "del")
  }

  # multinucleotide deletions 
  mbdel <- which((nchar(ref(x)) - nchar(alt(x))) > 1)
  mbdels <- length(mbdel)
  if (mbdels > 0) {
    mbstarts <- start(x)[mbdel] + 1
    mbends <- start(x)[mbdel] + (nchar(alt(x)) - nchar(ref(x)))[mbdel]
    nms[mbdel] <- paste0(mbstarts, "_", mbends, "del")
  }

  # get the inserted bases for a particular set of insertions:
  ins <- function(x, idx) psub(ref(x)[idx], alt(x)[idx])
  inses <- length(ins) 
  hgvsins <- function(x, idx) {
    paste0(start(x)[idx], "_", start(x)[idx]+1, "ins", ins(x, idx))
  }

  # indices of insertions (if any were detected)
  bpins <- which(nchar(alt(x)) > nchar(ref(x)))

  # rename insertions (length is not important) 
  if (length(bpins) > 0) nms[bpins] <- hgvsins(x, bpins)

  # if it's not an insertion or deletion, for now, it must be a substitution
  subs <- length(x) - (inses + sbdels + mbdels)
  
  # report 
  if (verbose) {
    if (subs > 0) message(plural(subs, "single-base subsitutions"))
    if (sbdels > 0) message(plural(sbdels, "single-base deletion"))
    if (mbdels > 0) message(plural(mbdels, "multi-base deletion"))
    if (inses > 0) message(plural(inses, "insertion"))
  }

  # adds prefix
  return(paste0(pre, nms))

}
