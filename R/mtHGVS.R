#' convert mitochondrial variant calls to HGVS format for naming
#'
#' @param mvr   an MVRanges
#'
#' @return      a properly named MVRanges
#'
#' @export
mtHGVS <- function(mvr) { 

  # baseline notation, assuming mostly SNVs 
  nms <- paste0(pos(mvr), ref(mvr), ">", alt(mvr))

  # deletions
  sbdel <- which((nchar(ref(mvr)) - nchar(alt(mvr))) == 1)
  mbdel <- which((nchar(ref(mvr)) - nchar(alt(mvr))) > 1)

  # single nucleotide deletions 
  if (length(sbdel) > 0) {
    nms[sbdel] <- paste0("g.", start(mvr)[sbdel] + 1, "del")
  }

  # multinucleotide deletions 
  if (length(mbdel) > 0) {
    mbstarts <- start(mvr)[mbdel] + 1
    mbends <- start(mvr)[mbdel] + (nchar(alt(mvr)) - nchar(ref(mvr)))[mbdel]
    nms[mbdel] <- paste0("g.", mbstarts, "_", mbends, "del")
  }

  # for insertion strings:
  psub <- function(x, y) { 
    stopifnot(length(x) == length(y))
    for (i in seq_along(x)) y[i] <- sub(paste0("^", x[i]), "", y[i])
    return(y)
  }

  # get the inserted bases for a particular set of insertions:
  ins <- function(mvr, idx) psub(ref(mvr)[idx], alt(mvr)[idx])
  hgvsins <- function(mvr, idx) {
    paste0("g.", start(mvr)[idx], "_", start(mvr)[idx]+1, "ins", ins(mvr, idx))
  }

  # indices of insertions (if any were detected)
  bpins <- which(nchar(alt(mvr)) > nchar(ref(mvr)))

  # rename insertions (length is not important) 
  if (length(bpins) > 0) nms[bpins] <- hgvsins(mvr, bpins)

  # done (for now)
  return(nms) 

}
