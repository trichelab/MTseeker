#' Mask haplogroup specific variants using haplogrep
#'
#' @name haploMask
#' 
#' @param mvr Input should be either an MVRanges or MVRangesList
#' @param fasta.output Path to a location to output the per-sample consensus variant fasta files
#' @param mask Whether to mask out the haplogroup specific SNPs
#' @param return.haplogroup Return the inferred haplogroup
#' @param ties Whether to use depth or weight to resolve variants to pass for haplogroup inference
#' @param override This will override the java check and proceed even if a java install isn't detected
#' @param java.path This is an option to specify a path to a java install (e.g. /usr/bin/java)
#'
#' @return An MVRanges or MVRangesList of masked haplogroup-specific variants
#'
#' @import utils
#' @import rtracklayer
#' @import GenomicRanges
#' 
#' @export
#'
#' @examples
#' 

haploMask <- function(mvr, fasta.output = NULL, mask = TRUE,
                      return.haplogroup = TRUE, ties = c("depth", "weight"),
                      override = FALSE, java.path = NULL) {

  #check if java is installed and stop if not
  if (is.null(java.path)) {
  java_install <- ifelse(length(system("java -version 2>&1", intern = TRUE)) == 3,
                         TRUE,
                         FALSE)
  } else {
    #user has provided a path to java install
    override <- TRUE
  }
  if (isFALSE(java_install) & isFALSE(override)) {
    stop("Couldn't find java. Make sure running 'java -version' doesn't return an error. Stopping.")
  }
  
  if (is(mvr, "MVRangesList")) {
    mvrl <- lapply(mvr, haploMask, fasta.output = fasta.output, mask = mask, return.haplogroup = return.haplogroup)
    return(suppressWarnings(MVRangesList(mvrl)))
    }
  
  #filter out non-PASSing variants
  #this should be done by default and not up to the user...
  #mvr <- filterMTvars(mvr)

  if (isEmpty(mvr)) return(mvr)
  
  #collect variants for haplogroup inference with haplogrep
  #this will fail with disjoint ranges
  #to deal intelligently, we need a whitelist of ranges/variants to pick from
  #otherwise we can just drop the range for haplogrep
  if (!isDisjoint(mvr)) {
    message("Found non-disjoint ranges...")
    #pull in whitelist
    data(haplomask_whitelist, package = "MTseeker")
    #find and check the non-disjoint ranges against the whitelist
    non_disjoint_ranges <- disjoin(mvr, with.revmap = TRUE)
    #subset the ranges to only keep the overlapping ranges
    disjoin_idxs <- lapply(mcols(non_disjoint_ranges)$revmap, function(x) {
      if (length(x) > 1) return(x)
      else {return(NA)}
    })
    disjoin_idxs <- unlist(disjoin_idxs[!is.na(disjoin_idxs)])
    #gather the disjoint ranges to query against whitelist
    disjoint_ranges <- mvr[disjoin_idxs,]
    
    message("Resolving variants for haplogroup inference...")
    #remove non-disjoint ranges from original object
    mvr <- mvr[-c(disjoin_idxs),]
    nd_snvs <- disjoint_ranges[isSNV(disjoint_ranges),]
    nd_ins <- disjoint_ranges[isInsertion(disjoint_ranges),]
    nd_del <- disjoint_ranges[isDeletion(disjoint_ranges),]
    #standard SNPs
    if (length(nd_snvs)) {
      ssnps <- subsetByOverlaps(nd_snvs, haplomask_whitelist$Standard_SNPs)
      #check for duplicate ranges
      if (length(ssnps)) {
        ssnps <- .resolveTies(ssnps, ties = ties)
        mvr <- sort(c(mvr, ssnps))
      }
    }
    #insertions
    #TODO: deal with .XC expansions to handle variable length insertions
    if (length(nd_ins)) {
      ins <- subsetByOverlaps(nd_ins, haplomask_whitelist$Insertions)
      #check for duplicate ranges
      if (length(ins)) {
        ins <- .resolveTies(ins, ties = ties)
        mvr <- sort(c(mvr, ins))
      }
    }
    #deletions
    if (length(nd_del)) {
      del <- subsetByOverlaps(nd_del, haplomask_whitelist$Deletions)
      #check for duplicate ranges
      if (length(del)) {
        del <- .resolveTies(del, ties = ties)
        mvr <- sort(c(mvr, del))
      }
    }
    #one last check for disjoint ranges
    #this can in fact arise if we have overlapping indels with high heteroplasmy
    #specifically this shows up in bulk samples
    #non-ideal way to deal is to just disjoin and warn the user
    #that some haplogroup variants still exist in the sample
    if (!isDisjoint(mvr)) {
      message("WARNING! WARNING!")
      message("Non-disjoint ranges still persist")
      message("Some haplogroup variants will not be masked")
      message("Likely due to high heteroplasmy in the sample")
      message("Disjoining and proceeding with arbitrary, variants to ensure it is disjoint")
      resolved.disjoin <- disjoin(mvr)
      mvr <- subsetByOverlaps(mvr, resolved.disjoin, type = "equal")
    }
  }
  consensus_fasta <- .injectVariantsIntoReference(mvr)
  
  #export the fasta files
  if (is.null(fasta.output)) fasta.output <- getwd()
  message("Writing fasta for ", unique(sampleNames(mvr)))
  fasta_output_loc <- paste0(fasta.output, "/", unique(sampleNames(mvr)), ".consensus.fasta")
  rtracklayer::export(consensus_fasta, 
                      fasta_output_loc)
  
  #run haplogrep with the extended report to return SNPs found to call haplogroup
  haplogrep_output_loc <- paste0(fasta.output, "/", unique(sampleNames(mvr)), ".consensus.haplogrep.txt")
  message("Running haplogrep on ", unique(sampleNames(mvr)))
  if (is.null(java.path)) {
    haplogrep_output <- system(paste0("java -jar ",
                                    system.file("extdata", "haplogrep-2.1.20.jar", package = "MTseeker"),
                                    " --in ", fasta_output_loc,
                                    " --format fasta",
                                    " --out ", haplogrep_output_loc,
                                    " --extend-report"))
  } else {
    haplogrep_output <- system(paste0(java.path, " -jar ",
                                      system.file("extdata", "haplogrep-2.1.20.jar", package = "MTseeker"),
                                      " --in ", fasta_output_loc,
                                      " --format fasta",
                                      " --out ", haplogrep_output_loc,
                                      " --extend-report"))
    }
  
  #read the output in and parse to mask
  message("Parsing haplogrep output for ", unique(sampleNames(mvr)))
  haplogrep_input <- read.delim(haplogrep_output_loc,
                                header = TRUE,
                                stringsAsFactors = FALSE)
  
  #return back the haplogroup
  if (return.haplogroup) {
    metadata(mvr)$haplogroup <- haplogrep_input$Haplogroup
    #add the haplogroup quality score
    metadata(mvr)$haplogroup.quality <- haplogrep_input$Quality
    #return the found polymorphisms
    metadata(mvr)$haplogroup.found.SNPs <- haplogrep_input$Found_Polys
    #return the input polymorphisms
    metadata(mvr)$haplogroup.input.SNPs <- haplogrep_input$Input_Sample
  }
  
  #check if no variants were found to mask from haplogroups
  if (is.na(haplogrep_input$Found_Polys)) return(mvr)
  
  #mask off the haplogroup-specific polymorphisms
  #match up input polymorphisms to found polymorphisms
  haplogroup_polys_to_mask <- strsplit(haplogrep_input$Found_Polys, split = " ")
  haplogroup_input_polys <- strsplit(haplogrep_input$Input_Sample, split = " ")
  #yuck... not sure why haplogrep does this... but this is critical
  haplogroup_polys_to_mask <- list(haplogroup_polys_to_mask[[1]][haplogroup_polys_to_mask[[1]] %in% haplogroup_input_polys[[1]]])
  haplogroup_polys_to_mask <- IRanges(start = as.numeric(sapply(haplogroup_polys_to_mask, function(x) gsub("([0-9]+).*$", "\\1", x))),
                                      end = as.numeric(sapply(haplogroup_polys_to_mask, function(x) gsub("([0-9]+).*$", "\\1", x))))
  haplogroup_polys_to_mask_alt <- strsplit(haplogrep_input$Found_Polys, split = " ")
  #yuck... not sure why haplogrep does this... but this is critical
  haplogroup_polys_to_mask_alt <- list(haplogroup_polys_to_mask_alt[[1]][haplogroup_polys_to_mask_alt[[1]] %in% haplogroup_input_polys[[1]]])
  haplogroup_polys_to_mask_alt <- sapply(haplogroup_polys_to_mask_alt, function(x) gsub("^([0-9]+)", "", x))
  
  #subset the existing mvr to get the ref seqs
  haplogroup_polys_to_mask.gr <- GRanges(seqnames = "chrM", ranges = haplogroup_polys_to_mask)
  haplogroup_polys_to_mask_ref <- subsetByOverlaps(mvr, haplogroup_polys_to_mask.gr)
  haplogroup_polys_to_mask_ref <- ref(haplogroup_polys_to_mask_ref[complete.cases(match(alt(haplogroup_polys_to_mask_ref), haplogroup_polys_to_mask_alt)),])
  snps_to_mask <- MVRanges(VRanges(seqnames = "chrM",
                                   ranges = haplogroup_polys_to_mask,
                                   ref = haplogroup_polys_to_mask_ref,
                                   alt = haplogroup_polys_to_mask_alt,
                                   totalDepth = 100,
                                   refDepth = 0,
                                   altDepth = 100,
                                   sampleNames = "Haplogrep"))
  if (mask) {
    masked <- subsetByOverlaps(mvr, snps_to_mask, invert = TRUE)
    return(masked)
    }
  masked <- snps_to_mask
  return(masked)
}

#helper function to resolve ties
.resolveTies <- function(mvr, ties = c("depth", "weight")) {
  #resolve each tie
  #currently depth is only supported
  #set up a way to filter
  metadata(mvr)$keep <- rep(NA, length(mvr))
  names(metadata(mvr)$keep) <- names(mvr)
  for (var in 1:length(mvr)) {
    if (ties == "depth") {
      #get equal ranges
      dup_ranges <- subsetByOverlaps(mvr,
                                     mvr[var,],
                                     type = "equal")
      #sort by the higher depth ranges
      if (length(unique(altDepth(dup_ranges))) > 1 & all(is.na(metadata(dup_ranges)$keep))) {
        higher_depth_var <- sort(altDepth(dup_ranges), decreasing = TRUE)
        metadata(mvr)$keep[dup_ranges@altDepth == higher_depth_var[1]] <- "keep"
      }
      
      if (length(unique(altDepth(dup_ranges))) == 1) {
        message("Couldn't resolve the tie based on depth.")
        message("Returning the first variant in the tie.")
        #this is a catch for iterating through the ranges to decompose them
        #ensures we don't set "keep" to both duplicate ranges and end up with non-disjoint ranges again...
        if (all(is.na(metadata(dup_ranges)$keep[names(dup_ranges)]))) {
          metadata(mvr)$keep[names(dup_ranges)[1]] <- "keep"
        }
      }
    }
  }
  resolved <- mvr[!is.na(metadata(mvr)$keep),]
  metadata(resolved) <- list()
  #one last check for disjoint ranges
  #this can in fact arise if we have overlapping indels with high heteroplasmy
  #specifically this shows up in bulk samples
  #non-ideal way to deal is to just disjoin and warn the user
  #that some haplogroup variants still exist in the sample
  if (!isDisjoint(resolved)) {
    message("WARNING! WARNING!")
    message("Non-disjoint ranges still persist")
    message("Some haplogroup variants will not be masked")
    message("Likely due to high heteroplasmy in the sample")
    message("Disjoining and proceeding with arbitrary, variants to ensure it is disjoint")
    resolved.disjoin <- disjoin(resolved)
    resolved <- subsetByOverlaps(resolved, resolved.disjoin, type = "equal")
  }
  return(resolved)
}

#bug in snpCall work around
.injectVariantsIntoReference <- function(x, ...) {
  supported <- c("rCRS")
  actual <- unique(genome(x))
  stopifnot(unique(genome(x)) %in% supported)
  data(rCRSeq, package="MTseeker")
  mvr <- x
  alts <- DNAStringSet(replaceAt(rCRSeq[[1]], ranges(mvr), alt(mvr)))
  names(alts) <- actual # genome 
  return(alts)
}
