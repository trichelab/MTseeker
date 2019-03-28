#' Mask haplogroup specific variants using haplogrep
#'
#' @name haploMask
#' 
#' @param mvr Input should be either an MVRanges or MVRangesList
#' @param fasta.output Path to a location to output the per-sample consensus variant fasta files
#' @param mask Whether to mask out the haplogroup specific SNPs
#' @param return.haplogroup Return the inferred haplogroup
#' @param override This will override the java check and proceed even if a java install isn't detected
#'
#' @return An MVRanges or MVRangesList of masked haplogroup-specific variants
#' 
#' @import rtracklayer
#' @import GenomicRanges
#' 
#' @export
#'
#' @examples
#' 

haploMask <- function(mvr, fasta.output = NULL, mask = TRUE, return.haplogroup = TRUE, override = FALSE) {

  #check if java is installed and stop if not
  java_install <- ifelse(length(system("java -version 2>&1", intern = TRUE)) == 3,
                         TRUE,
                         FALSE)
  if (isFALSE(java_install) & isFALSE(override)) {
    stop("Couldn't find java. Make sure running 'java -version' doesn't return an error. Stopping.")
  }
  
  if (is(mvr, "MVRangesList")) {
    mvrl <- lapply(mvr, haploMask, fasta.output = fasta.output, mask = mask, return.haplogroup = return.haplogroup)
    return(suppressWarnings(MVRangesList(mvrl)))
    }
  
  #collect variants for haplogroup inference with haplogrep
  consensus_fasta <- consensusString(mvr)
  
  #export the fasta files
  if (is.null(fasta.output)) fasta.output <- getwd()
  message("Writing fasta for ", unique(sampleNames(mvr)))
  fasta_output_loc <- paste0(fasta.output, "/", unique(sampleNames(mvr)), ".consensus.fasta")
  rtracklayer::export(consensus_fasta, 
                      fasta_output_loc)
  
  #run haplogrep with the extended report to return SNPs found to call haplogroup
  haplogrep_output_loc <- paste0(fasta.output, "/", unique(sampleNames(mvr)), ".consensus.haplogrep.txt")
  message("Running haplogrep on ", unique(sampleNames(mvr)))
  haplogrep_output <- system(paste0("java -jar ",
                                    system.file("extdata", "haplogrep-2.1.20.jar", package = "MTseeker"),
                                    " --in ", fasta_output_loc,
                                    " --format fasta",
                                    " --out ", haplogrep_output_loc,
                                    " --extend-report"))
  
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
  }
  
  #check if no variants were found to mask from haplogroups
  if (is.na(haplogrep_input$Found_Polys)) return(mvr)
  
  #mask off the haplogroup-specific polymorphisms
  haplogroup_polys_to_mask <- strsplit(haplogrep_input$Found_Polys, split = " ")
  haplogroup_polys_to_mask <- IRanges(start = as.numeric(sapply(haplogroup_polys_to_mask, function(x) gsub("([0-9]+).*$", "\\1", x))),
                                      end = as.numeric(sapply(haplogroup_polys_to_mask, function(x) gsub("([0-9]+).*$", "\\1", x))))
  haplogroup_polys_to_mask_alt <- strsplit(haplogrep_input$Found_Polys, split = " ")
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