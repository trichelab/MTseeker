#' Annotate AA changes in MT variants
#'
#' @name getProteinImpact
#'
#' @param mvr         An MVRangesList or MVRanges object
#' @param coding      TRUE implies look at only coding regions
#' @param parallel    Whether to run things in parallel
#' 
#'
#' @return            Annotated variants
#'
#' @import            jsonlite
#' @examples
#' 
#' 
#' @export

getProteinImpact <- function(mvr, coding=TRUE, parallel=FALSE, cores=1) {

  if (genome(mvr) != "rCRS") {
    mvr <- decomposeAndCalcConsequences(mvr, coding=coding)
    
    if (is(mvr, "MVRanges")) mvr$apogee <- FALSE
    else {
      for (i in 1:length(mvr)) {
        mcols(mvr[[i]])$apogee <- FALSE
      }
    }
    
    # When using decompAndCalcCons, if there are no AA changes, it will set it to ""
    # Setting it to NA if someone wants to remove synonymous mutation later on
    if (length(which(mvr$AAchange == "")) > 0) mvr[which(mvr$AAchange == "")]$AAchange <- NA_character_
    if (length(which(!mvr$apogee)) > 0)  mvr[which(!mvr$apogee)]$impacted.gene <- paste0("MT-", mvr[which(!mvr$apogee)]$impacted.gene) 
    
    # Adds a column determing type of mutation
    # Missense, synonymous, nonsense, insertion, deletion, or frameshift
    mvr <- .typeMutation(mvr)
    
    return(mvr) 
  }
  
  # Set the number of cores if running in parallel
  if (parallel) {
    message("Parallel execution is currently disabled")
  }
  if (is(mvr, "MVRangesList")) {
    MVRangesList(lapply(mvr, getProteinImpact, coding=coding))
  }
  
  stopifnot(is(mvr, "MVRanges"))
  
  message("Processing consequences for ", sampleNames(mvr)@values)

  # Right now can only handle variants in coding regions
  if (coding) mvr <- .getCodingRegions(mvr)
  
  # Return empty MVRanges if there are none in coding regions
  if (length(mvr) == 0) return(mvr)
  
  # Output that we want to have by the end of this
  mvr$AAchange <- NA_character_
  mvr$impacted.gene <- NA_character_
  mvr$overlapGene <- NA_character_
  mvr$apogee <- TRUE

  for (i in 1:length(mvr)) {

    # Get the url and scrape the information from mitimpact
    pos <- paste0(start(mvr[i]),  "-", end(mvr[i]))
    url <- paste0(paste("http://mitimpact.css-mendel.it", "api", "v2.0", 
                        "genomic_position/", sep="/"), pos)
    res <- as.data.frame(read_json(url, simplifyVector=TRUE)$variants, stringsAsFactors=F)
    
    # If mitimpact gives a valid result
    # Mitimpact does not handle insertions and deletion?
    if (nrow(res) > 0 && !(grepl("ins|del", names(mvr[i])))) {
      
      # Rename the variant and paste together a simpler version of the consequences
      res$genomic <- with(res, paste0("m.", Start, Ref, ">", Alt))
      res$consequence <- with(res, paste0(AA_ref,AA_position,AA_alt))
      
      # The information we are interested in
      columns <- c("genomic","protein","Start", "Ref","Alt","Codon_substitution","consequence",
                   "Mitomap_Phenotype","Mitomap_Status","Gene_symbol",
                   "OXPHOS_complex","Consequence","APOGEE_boost_consensus")
      columnsPresent <- intersect(columns, names(res))
      mitOut <- res[, columnsPresent]
      
      # Sometimes there are multiple variants that can be found at one position
      # Find the one that aligns with the variant we are looking at
      index <- which(mitOut$Alt == alt(mvr[i]))
      
      # If none of the results from apogee have the same alt and ref sequence
      # Go back to the old way of doing things
      if (length(index) == 0) {
        mvr[i] <- decomposeAndCalcConsequences(mvr[i], coding=coding)
        mvr[i]$apogee <- FALSE
      } 
      
      # If apogee gives a result that we want
      else {
        
        if (mitOut$Ref[index] != ref(mvr[i])) {
          message(paste0("Reference disparity between APOGEE and pileup for variant: ", names(mvr[i]), " in sample ", mvr[i]$sampleNames))
        }
        
        out <- mitOut[index,]
        mvr[i]$AAchange <- out$consequence
        mvr[i]$impacted.gene <- out$Gene_symbol
      }
    }
    
    # If mitimpact does not give a valid result, go back to the old way
    else {
      mvr[i] <- decomposeAndCalcConsequences(mvr[i], coding=coding)
      mvr[i]$apogee <- FALSE
    }
    
  } #for
 
  # When using decompAndCalcCons, if there are no AA changes, it will set it to ""
  # Setting it to NA if someone wants to remove synonymous mutation later on
  if (length(which(mvr$AAchange == "")) > 0) mvr[which(mvr$AAchange == "")]$AAchange <- NA_character_
  if (length(which(!mvr$apogee)) > 0)  mvr[which(!mvr$apogee)]$impacted.gene <- paste0("MT-", mvr[which(!mvr$apogee)]$impacted.gene) 
  
  # Adds a column determing type of mutation
  # Missense, synonymous, nonsense, insertion, deletion, or frameshift
  mvr <- .typeMutation(mvr)

  return(mvr) 
}

# helper function to subset ranges to just coding space
.getCodingRegions <- function(mvr) {
  
  # rCRS only, for the time being 
  stopifnot(unique(genome(mvr)) == "rCRS")
  
  # get mtGenes if needed 
  gr <- genes(mvr)
  stopifnot(unique(genome(gr)) == "rCRS")
  
  mvr <- MVRanges(subsetByOverlaps(mvr, gr, type="within"), coverage=genomeCoverage(mvr))
  
  # subset the variants to those that overlap the target GRanges and are canon
  if (length(mvr)) {
    
    #drop anything that has an N base.. this also looks like a weird bug?
    mvr <- mvr[!grepl("N", mvr@alt),]
    
    #check again whether we've now cleared out all the variants
    #return an empty ranges if we have
    if (length(mvr) == 0) {
      mvr <- MVRanges(subsetByOverlaps(mvr, gr, type="within"))
    }
  } 
  return(mvr)
}

# helper fn to determine type of mutation
.typeMutation <- function(mvr) {
  
  mvr$typeMut <- NA_character_

  # Determine the overaching consequences (type of mutation)
  # SNP
  snp <- grep(">", names(mvr))
  if (length(snp) != 0) {
    
    # Synonymous : no AA change overall
    syn <- which(is.na(mvr[snp]$AAchange))
    mvr[snp]$typeMut[syn] <- "synonymous"
    
    # Missense or nonsense
    missOrNon <- which(!is.na(mvr[snp]$AAchange))
    for (i in 1:length(missOrNon)) {
      index <- missOrNon[i]
      
      # Want to split up the consequences ex. V12I to becomes "V" and "I" 
      # Where V is the ref AA and I is the alt AA to determine what change took place
      aaChange <- unlist(strsplit(gsub("\\d+", " ", mvr[snp]$AAchange[index]), " "))
      
      # If the * was already a part of the reference sequence
      # Then it is not a nonsense mutation
      if ( ("*" %in% aaChange[2]) && !("*" %in% aaChange[1]) ) mvr[snp]$typeMut[index] <- "nonsense"
      else mvr[snp]$typeMut[index] <- "missense"
    }
    
  }
  
  ins <- grep("ins", names(mvr))
  if (length(ins) != 0) {
    
    insLength <- nchar(alt(mvr[ins])) - 1
    
    mvr[ins]$typeMut[which(insLength %% 3 == 0)] <- "insertion"
    mvr[ins]$typeMut[which(!(insLength %% 3 == 0))] <- "frameshift"
    
  }
  
  del <- grep("del", names(mvr))
  if (length(del) != 0) {
    
    delLength <- nchar(ref(mvr[del])) - 1
    
    mvr[del]$typeMut[which(delLength %% 3 == 0)] <- "deletion"
    mvr[del]$typeMut[which(!(delLength %% 3 == 0))] <- "frameshift"
  }
  
  return(mvr)
}
