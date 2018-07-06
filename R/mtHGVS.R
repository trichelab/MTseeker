#' convert mitochondrial variant calls to HGVS format for naming
#'
#' @param x     an MVRanges (or, in a pinch, a GRanges)
#'
#' @return      a properly named MVRanges (or, in a pinch, GRanges)
#'
#' @export
mtHGVS <- function(x) { 

  # for plotting purposes
  if (is(x, "GRanges")) {
    return(paste0("m.", sapply(apply(cbind(start(x), end(x)), 1, unique), 
                               paste, collapse="-")))
  }

  # baseline notation, assuming mostly SNVs 
  nms <- paste0(gsub("\\-", "_", pos(x)), ref(x), ">", alt(x))

  # deletions
  sbdel <- which((nchar(ref(x)) - nchar(alt(x))) == 1)
  mbdel <- which((nchar(ref(x)) - nchar(alt(x))) > 1)

  # single nucleotide deletions 
  if (length(sbdel) > 0) {
    nms[sbdel] <- paste0(start(x)[sbdel] + 1, "del")
  }

  # multinucleotide deletions 
  if (length(mbdel) > 0) {
    mbstarts <- start(x)[mbdel] + 1
    mbends <- start(x)[mbdel] + (nchar(alt(x)) - nchar(ref(x)))[mbdel]
    nms[mbdel] <- paste0(mbstarts, "_", mbends, "del")
  }

  # for insertion strings:
  psub <- function(x, y) { 
    stopifnot(length(x) == length(y))
    for (i in seq_along(x)) y[i] <- sub(paste0("^", x[i]), "", y[i])
    return(y)
  }

  # get the inserted bases for a particular set of insertions:
  ins <- function(x, idx) psub(ref(x)[idx], alt(x)[idx])
  hgvsins <- function(x, idx) {
    paste0(start(x)[idx], "_", start(x)[idx]+1, "ins", ins(x, idx))
  }

  # indices of insertions (if any were detected)
  bpins <- which(nchar(alt(x)) > nchar(ref(x)))

  # rename insertions (length is not important) 
  if (length(bpins) > 0) nms[bpins] <- hgvsins(x, bpins)

  # prepend m.
  return(paste0("m.", nms))

}
