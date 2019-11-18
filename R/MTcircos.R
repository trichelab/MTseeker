#' plot a canonical human (or, in principle, any) mitochondrial genome 
#'
#' The default font sizes, orientations, etc. are optimized for a "cold" start;
#' if you want to fiddle with the details, crack open the code and modify it...
#' or alternatively, add sectors/dendrograms inside of this "framed" version.
#'
#' 
#' @param variants  optional MVRanges or MVRangesList to split by strand & plot
#' @param outside   optional MVRanges or MVRangesList to plot outside the circle
#' @param inside    optional MVRanges or MVRangesList to plot inside the circle
#' @param outcol    optional color assignment function or matrix for outside
#' @param incol     optional color assignment function or matrix for inside
#' @param anno      a GRanges (optional, defaults to mtAnno.rCRS if none given)
#' @param how       optional specification for how to plot multiple samples
#' @param ...       other arguments to pass on to called functions
#' 
#' @return          invisibly, a list: `anno` (data.frame) + `pfun` (panel.fun)
#'
#' @import VariantTools
#' @import rtracklayer
#' @import circlize
#' @import viridis
#'
#' @importFrom grDevices col2rgb rgb adjustcolor
#' @importFrom graphics legend
#'
#' @aliases genMTcircos initMTcircos genesMTcircos
#'
#' @examples 
#'
#' library(MTseekerData)
#' BAMdir <- system.file("extdata", "BAMs", package="MTseekerData")
#' BAMs <- paste0(BAMdir, "/", list.files(BAMdir, pattern="^pt.*bam$"))
#' (mvr <- filterMTvars(pileupMT(BAMs[1], ref="rCRS")))
#' MTcircos(mvr)
#' title("Mitochondrial variants in a leukemia patient's blast cell")
#' 
#' @export 
MTcircos <- function(variants=NULL, outside=NULL, inside=NULL, outcol=NULL, 
                     incol=NULL, anno=NULL, how=c("matrix", "AA"), ...) {
  circos.clear() 
  
  if (length(how) > 1) {
    how <- "matrix"
  }
  
  if (how == "AA") {
    
    if (is(variants, "MVRanges")) {
      if (!("typeMut" %in% names(mcols(variants)))) {
        message("Must run getProteinImpact() before plottting AA changes")
        stop()
      }
    }
    

    else if (is(variants, "MVRangesList")) {
      if(!"typeMut" %in% names(mcols(variants[[1]]))) {
        message("Must run getProteinImpact() before plottting AA changes")
        stop()
      }
    }
    
  }
  
  anno <- initMTcircos(variants)
  dat <- data.frame(name=names(anno), start=start(anno), end=end(anno))
  
  if (!is.null(variants) && length(variants) != 0) {
    message("Splitting variants by strand...")
    stranded <- byStrand(variants, anno)
    message("Replacing `outside` with heavy-strand variants...")
    outside <- stranded$heavy
    message("Replacing `inside` with light-strand variants...")
    inside <- stranded$light
  } 

  # Outside track
  if (!isEmpty(outside)) {
    
    # Color code according to AA changes
    if (how == "AA") {
      bed1 <- .AAmakeColoredMatrix(outside)
      vafBed1 <- .vafMatrix(bed1, outside, how)
    }
    
    else {
      # Color coding for the variants
      # del=blue, SNV=black, ins=red
      bed1 <- .makeColoredMatrix(outside)
      vafBed1 <- .vafMatrix(bed1, outside, how)
    }
    
    circos.genomicHeatmap(bed1, outcol, line_col=.colorCode(bed1$chr), 
                          col = vafBed1, track.margin=c(0,0), 
                          side="outside", border=NA,
                          line_lwd=2) 
  }
  else {
    circos.track(track.height=0.15, ylim=c(0,1), bg.border=NA)
  }
  
  # main track, gene names and such
  res <- genesMTcircos(variants, anno, legends=T)

  # Inside track
  if (!isEmpty(inside)) {
    
    # Color code according to AA changes
    if (how == "AA") {
      bed2 <- .AAmakeColoredMatrix(inside)
      vafBed2 <- .vafMatrix(bed2, inside, how)
    }
    
    else {
      # del=blue, SNV=black, ins=red
      bed2 <- .makeColoredMatrix(inside)
      vafBed2 <- .vafMatrix(bed2, inside, how)
    }
    
    circos.genomicHeatmap(bed2, outcol, line_col=.colorCode(bed2$chr), col = vafBed2,
                          track.margin=c(0,0), side="inside", border=NA,
                          line_lwd=2) 
  }
  else {
    circos.track(track.height=0.15, ylim=c(0,1), bg.border=NA)
  }
  
  # Color code according to AA changes
  if (how == "AA") {
    legend("topright", title="AA Change", ncol=2,
           legend=c("Missense", "Nonsense", "Synonymous", "Insertion", "Deletion", "Frameshift"), 
           col=viridis(6), pch=15, cex=0.6)
    
  }

  
  else {
    legend("topright", title="Variant Type",
           legend=c("Insertion", "Deletion", "SNV"), col=c("red", "blue", "black"), 
           pch=15, cex=0.8)
  }
 
  
  invisible(res)
  
}


# helper fn
.makeColoredMatrix <- function(mvr) {

  # Making a colored matrix
  if (is(mvr, "MVRangesList")) {
    allNames <- lapply(mvr, names)
    allNames <- unlist(unname(allNames))
    numSamples <- length(mvr)
  } else {
    allNames <- make.unique(names(mvr))
    numSamples <- 1
  }

  rowNam <- c("chr", "start", "end")

  if (numSamples == 1) {
    rowNam <- append(rowNam, unique(as.character(sampleNames(mvr))))
  } else {
    rowNam <- append(rowNam, names(mvr))
  }
  
  m <- matrix(0, ncol = length(rowNam), nrow = length(allNames))
  typeDF <- data.frame(m)

  names(typeDF) <- rowNam
  rownames(typeDF) <- make.unique(allNames)

  # Figure out which of the variants in a sample is SNV, ins, del
  # Assign color according to the unique values set for each type of variant
  if (is(mvr, "MVRanges")) {
    
    ins <- grep("ins", names(mvr))
    del <- grep("del", names(mvr))
    snv <- grep(">", names(mvr))

    typeDF[,4][ins] <- 1
    typeDF[,4][snv] <- 2
    typeDF[,4][del] <- 3

  } else {
    for (i in seq_len(numSamples)) {
      
      #if(length(mvr[[i]]) == 0) next
      
      rowInd <- which(allNames %in% names(mvr[[i]]))
      rowOverlapNames <- allNames[rowInd]
      
      snvs <- rowOverlapNames[grep(">", rowOverlapNames)]
      ins <- rowOverlapNames[grep("ins", rowOverlapNames)]
      dels <- rowOverlapNames[grep("del", rowOverlapNames)]
      
      if (length(ins) > 0) typeDF[ins,][,i + 3] <- 1
      if (length(snvs) > 0) typeDF[snvs,][,i + 3] <- 2
      if (length(dels) > 0) typeDF[dels,][,i + 3] <- 3
      # del=blue, SNV=black, ins=red

    }
  }
  
  # Get the start and end positions
  if (is(mvr, "MVRanges")) {
    typeDF$start <- start(mvr)
    typeDF$end <- end(mvr)
  }
  
  else {
    # Get the start and end position for each variants
    for (j in seq_len(numSamples)) {
      
      #if(length(mvr[[j]]) == 0) next
      
      hits <- which(typeDF[,j + 3] != 0)
      varNames <- allNames[hits]
      
      vars <- mvr[[j]][varNames]
      varStart <- start(vars)
      varEnd <- end(vars)
      
      typeDF$start[hits] <- varStart
      typeDF$end[hits] <- varEnd
    }
  }

  # Get the names of the genes each variant is found in
  anno <- suppressMessages(getAnnotations(mvr))
  ov <- findOverlaps(IRanges(typeDF$start, typeDF$end), ranges(anno))
  if (length(ov) > 0) { 
    # If there are overlapping genes
    # Only list the first gene it overlaps in
    firstOv <- ov[match(unique(queryHits(ov)), queryHits(ov)), ]
    typeDF$chr <- names(anno)[subjectHits(firstOv)]
  }

  return(typeDF)
}

.vafMatrix <- function(bed, mvr, how) {

  vafs <- bed
  vafs[,4:ncol(vafs)] <- 0


  # Get the VAF for each nonzero element of the matrix
  if (is(mvr, "MVRanges")) {
    vafs[,4] <- mvr$VAF
  } else {
    # Get the VAF for each variants
    for (j in seq_len(ncol(bed) - 3)) {
      
      hits <- which(bed[,j + 3] != 0)
      varNames <- row.names(bed)[hits]
      
      vars <- mvr[[j]][varNames]
      varVAF <- vars$VAF
      
      vafs[,j + 3][hits] <- varVAF
    }
  }


  # Nonzero elements
  nonzero <- which(bed[,4:ncol(bed)] !=0, arr.ind=T)
  
  # Color according to AA change
  if (how == "AA") {
    
    cols <- viridis(6)
    
    if (is(mvr, "MVRanges")) {
      
      bed[,4][which(bed[,4] == 1)] <- cols[1]
      bed[,4][which(bed[,4] == 2)] <- cols[2]
      bed[,4][which(bed[,4] == 3)] <- cols[3]
      bed[,4][which(bed[,4] == 4)] <- cols[4]
      bed[,4][which(bed[,4] == 5)] <- cols[5]
      bed[,4][which(bed[,4] == 6)] <- cols[6]
      
      for (i in seq_along(mvr)) {
        bed[,4][i] <- adjustcolor(col = bed[,4][i],alpha.f = vafs[,4][i])
      }
    } else {

      for (i in seq_len(ncol(bed) - 3)) {
        
        bed[,i + 3][which(bed[,i + 3] == 1)] <- cols[1]
        bed[,i + 3][which(bed[,i + 3] == 2)] <- cols[2]
        bed[,i + 3][which(bed[,i + 3] == 3)] <- cols[3]
        bed[,i + 3][which(bed[,i + 3] == 4)] <- cols[4]
        bed[,i + 3][which(bed[,i + 3] == 5)] <- cols[5]
        bed[,i + 3][which(bed[,i + 3] == 6)] <- cols[6]
        
      }

      # Add transparency
      for (k in seq_len(nrow(nonzero))) {
        
        # 1 is the row
        # 2 is the column
        rows <- nonzero[k,][1]
        cols <- nonzero[k,][2]
        
        bed[rows,][cols + 3] <- adjustcolor(col = bed[rows,][cols + 3], 
                                            alpha.f = vafs[rows,][cols + 3])
      }
      
    }
    
  } else {
    
    # Only color code snvs, ins, dels
    if (is(mvr, "MVRanges")) {

      bed[,4][which(bed[,4] == 1)] <- "red"
      bed[,4][which(bed[,4] == 2)] <- "black"
      bed[,4][which(bed[,4] == 3)] <- "blue"
      
      for (i in seq_along(mvr)) {
        bed[,4][i] <- adjustcolor(col = bed[,4][i],alpha.f = vafs[,4][i])
      }
    } else {
      for (k in seq_len(nrow(nonzero))) {
        
        # 1 is the row
        # 2 is the column
        rows <- nonzero[k,][1]
        cols <- nonzero[k,][2]
        
        #col_fun = colorRamp2(c(0, 1, 2, 3), c("white", "red", "black", "blue"))
        if (bed[rows,][cols + 3] == 1) bed[rows,][cols + 3] <- "red"
        else if (bed[rows,][cols + 3] == 2) bed[rows,][cols + 3] <- "black"
        else if (bed[rows,][cols + 3] == 3) bed[rows,][cols + 3] <- "blue"
        
        bed[rows,][cols + 3] <- adjustcolor(col = bed[rows,][cols + 3], 
                                            alpha.f = vafs[rows,][cols + 3])
      }
    }

  }

  bed <- bed[, setdiff(seq_len(ncol(bed)), c(1,2,3))]
  if (is(mvr, "MVRanges")) bed <- as.matrix(bed)
  return(bed)
}

.AAmakeColoredMatrix <- function(mvr) {
  
  # Making a colored matrix
  if (is(mvr, "MVRangesList")) {
    allNames <- lapply(mvr, names)
    allNames <- make.unique(unlist(unname(allNames)))
    numSamples <- length(mvr)
  } else {
    allNames <- make.unique(names(mvr))
    numSamples <- 1
  }
  
  rowNam <- c("chr", "start", "end")
  
  if (numSamples == 1) rowNam <- append(rowNam, unique(as.character(sampleNames(mvr))))
  else rowNam <- append(rowNam, names(mvr))
  
  m <- matrix(0, ncol = length(rowNam), nrow = length(allNames))
  typeDF <- data.frame(m)
  
  names(typeDF) <- rowNam
  rownames(typeDF) <- make.unique(allNames)
  
  
  # Assign values to the type of mutations
  # nonsense, synonymous, insertion, deletion, etc.
  if (is(mvr, "MVRanges")) {
    
    missense <- grep("missense", mvr$typeMut)
    nonsense <- grep("nonsense", mvr$typeMut)
    synonymous <- grep("synonymous", mvr$typeMut)
    insertion <- grep("insertion", mvr$typeMut)
    deletion <- grep("deletion", mvr$typeMut)
    frameshift <- grep("frameshift", mvr$typeMut)
    
    typeDF[,4][missense] <- 1
    typeDF[,4][nonsense] <- 2
    typeDF[,4][synonymous] <- 3
    typeDF[,4][insertion] <- 4
    typeDF[,4][deletion] <- 5
    typeDF[,4][frameshift] <- 6
    
  }
  
  else {
    
    # Figure out which of the variants in a sample is which type of variants
    # Assign color according to the unique values set for each type of variant
    for (i in seq_len(numSamples)) {

      rowInd <- which(allNames %in% names(mvr[[i]]))
      rowOverlapNames <- allNames[rowInd]

      missense <- rowInd[grep("missense", mvr[[i]][rowOverlapNames,]$typeMut)]
      nonsense <- rowInd[grep("nonsense", mvr[[i]][rowOverlapNames,]$typeMut)]
      synonymous <- rowInd[grep("synonymous", mvr[[i]][rowOverlapNames,]$typeMut)]
      insertion <- rowInd[grep("insertion", mvr[[i]][rowOverlapNames,]$typeMut)]
      deletion <- rowInd[grep("deletion", mvr[[i]][rowOverlapNames,]$typeMut)]
      frameshift <- rowInd[grep("frameshift", mvr[[i]][rowOverlapNames,]$typeMut)]
      
      typeDF[,i+3][missense] <- 1
      typeDF[,i+3][nonsense] <- 2
      typeDF[,i+3][synonymous] <- 3
      typeDF[,i+3][insertion] <- 4
      typeDF[,i+3][deletion] <- 5
      typeDF[,i+3][frameshift] <- 6
    
    }
  }
  
  # Get the start and end positions for each variant
  if (is(mvr, "MVRanges")) {

    typeDF$start <- start(mvr)
    typeDF$end <- end(mvr)
    
  }
  
  else {
    
    # Get the start and end position for each variants
    for (j in seq_len(numSamples)) {
      
      hits <- which(typeDF[,j + 3] != 0)
      varNames <- allNames[hits]
      
      vars <- mvr[[j]][varNames]
      varStart <- start(vars)
      varEnd <- end(vars)
      
      typeDF$start[hits] <- varStart
      typeDF$end[hits] <- varEnd
    }
    
  }
  
  
  # Get the names of the genes each variant is found in
  anno <- suppressMessages(getAnnotations(mvr))
  ov <- findOverlaps(IRanges(typeDF$start, typeDF$end), ranges(anno))
  
  # If there are overlapping genes
  # Only list the first gene it overlaps in
  firstOv <- ov[match(unique(queryHits(ov)), queryHits(ov)), ]
  
  typeDF$chr <- names(anno)[subjectHits(firstOv)]
  
  return(typeDF)
}
