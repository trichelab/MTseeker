#' List all defined methods for an S4 class (or classes, if you must)
#' 
#' Mostly for debugging and making architectural choices, e.g. about coverage()
#' 
#' Note: this is borrowed from Hadley, who borrowed it from a BioC workshop!
#' 
#' @param ...     name[s] of class[es] (please note, results will be union'ed)
#' 
#' @return        methods for the class[es], union'ed into a character vector
#' 
#' @import        methods
#' 
#' @examples
#' s4methods("MVRangesList")
#' s4methods("MAlignmentsList")
#'
#' @export
s4Methods <- function(...) {
  methods <- showMethods(classes=as.character(...), printTo=FALSE)
  methods <- methods[grep("^Function:", methods)]
  sapply(strsplit(methods, " "), "[", 2)
}
