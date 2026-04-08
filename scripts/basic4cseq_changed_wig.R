printWigFile <- function(expData, wigFileName = "output.wig", fixedSpan = TRUE, headerUCSC = "", useOnlyIndex = FALSE) {
  
  fragmentData = rawFragments(expData)
  chromosomes = unique(as.vector(fragmentData$chromosomeName))
  
  # prepare chromosome names for wig --> even if "1"..."Y" have been used for the library or in the reference genome, 
  # one may wish to print the wig with "chr1"..."chrY" notations
  chromosomeNamesWig = gsub("chr", "", as.character(chromosomes))
  if (!useOnlyIndex) {
    chromosomeNamesWig = paste("chr", chromosomeNamesWig, sep = "")
  }
  
  wigFileContent = NULL
  
  if (headerUCSC != "") {
    wigFileContent = headerUCSC
  }
  
  for (i in 1:length(chromosomes)) {
    
    fragsChromosome = subset(fragmentData, fragmentData$chromosomeName == chromosomes[i])
    chromosomeWig = giveWigDataChromosome(fragsChromosome, readLength(expData), chromosomeNamesWig[i], fixedSpan)
    
    wigFileContent = c(wigFileContent, chromosomeWig)
  }
  
  fileConn <- file(wigFileName)
  
  writeLines(wigFileContent, fileConn)
  
  close(fileConn)
}

giveWigDataChromosome <- function(fragmentDataChromosome, readLength, chromosomeID, fixedSpan = TRUE) {
  
  if (!fixedSpan) {
    # variable step
    wigHeader = paste("variableStep chrom=", chromosomeID, sep = "")
  } else {
    # variable step with fixed span (read length)
    wigHeader = paste("variableStep chrom=", chromosomeID, " span=", readLength, sep = "")
  }
  
  # get reads on valid unique fragments
  fragmentDataLeftValid = subset(fragmentDataChromosome, fragmentDataChromosome$leftFragEndValid == TRUE)
  fragmentDataRightValid = subset(fragmentDataChromosome, fragmentDataChromosome$rightFragEndValid == TRUE)
  
  fragsLeft = cbind(fragmentDataLeftValid$fragmentStart, fragmentDataLeftValid$leftFragEndReads)
  fragsRight = cbind(fragmentDataRightValid$fragmentEnd - readLength + 1, fragmentDataRightValid$rightFragEndReads)
  
  if (!fixedSpan) {
    fragsLeftPos = NULL
    fragsLeftReads = NULL
    fragsRightPos = NULL
    fragsRightReads = NULL
    for (i in 1:readLength) {
      fragsLeftPos = c(fragsLeftPos, fragsLeft[,1]+i-1)
      fragsLeftReads = c(fragsLeftReads, fragsLeft[,2])
      fragsRightPos = c(fragsRightPos, fragsRight[,1]+i-1)
      fragsRightReads = c(fragsRightReads, fragsRight[,2])
    }
    fragsLeft = cbind(fragsLeftPos, fragsLeftReads)
    fragsRight = cbind(fragsRightPos, fragsRightReads)
  }
  
  fragEnds = rbind(fragsLeft, fragsRight)
  fragEnds = fragEnds[order(fragEnds[,1]),]
  
  wigData = c(wigHeader, paste(fragEnds[,1], fragEnds[,2], sep = " "))
  
  return(wigData)
}




setMethod("printWigFile",
          signature=signature(expData="Data4Cseq"),
          printWigFile)


setMethod("giveWigDataChromosome",
          signature=signature(fragmentDataChromosome="data.frame", readLength="numeric", chromosomeID="character"),
          giveWigDataChromosome)

printWigFile2 <- function(expData, wigFileName = "output.wig", fixedSpan = TRUE, headerUCSC = "", useOnlyIndex = FALSE) {
  
  fragmentData = rawFragments(expData)
  chromosomes = unique(as.vector(fragmentData$chromosomeName))
  
  # prepare chromosome names for wig --> even if "1"..."Y" have been used for the library or in the reference genome, 
  # one may wish to print the wig with "chr1"..."chrY" notations
  chromosomeNamesWig = gsub("chr", "", as.character(chromosomes))
  if (!useOnlyIndex) {
    chromosomeNamesWig = paste("chr", chromosomeNamesWig, sep = "")
  }
  
  wigFileContent = NULL
  
  if (headerUCSC != "") {
    wigFileContent = headerUCSC
  }
  
  for (i in 1:length(chromosomes)) {
    
    fragsChromosome = subset(fragmentData, fragmentData$chromosomeName == chromosomes[i])
    chromosomeWig = giveWigDataChromosome2(fragsChromosome, readLength(expData), chromosomeNamesWig[i], fixedSpan)
    
    wigFileContent = c(wigFileContent, chromosomeWig)
  }
  
  fileConn <- file(wigFileName)
  
  writeLines(wigFileContent, fileConn)
  
  close(fileConn)
}

giveWigDataChromosome2 <- function(fragmentDataChromosome, readLength, chromosomeID, fixedSpan = TRUE) {
  
  if (!fixedSpan) {
    # variable step
    wigHeader = paste("variableStep chrom=", chromosomeID, sep = "")
  } else {
    # variable step with fixed span (read length)
    wigHeader = paste("variableStep chrom=", chromosomeID, " span=", readLength, sep = "")
  }
  
  # get reads on valid unique fragments
  #fragmentDataLeftValid = subset(fragmentDataChromosome, fragmentDataChromosome$leftFragEndValid == TRUE)
  #fragmentDataRightValid = subset(fragmentDataChromosome, fragmentDataChromosome$rightFragEndValid == TRUE)
  
  fragsLeft = cbind(fragmentDataChromosome$fragmentStart, fragmentDataChromosome$leftFragEndReads)
  fragsRight = cbind(fragmentDataChromosome$fragmentEnd - readLength + 1, fragmentDataChromosome$rightFragEndReads)
  
  if (!fixedSpan) {
    fragsLeftPos = NULL
    fragsLeftReads = NULL
    fragsRightPos = NULL
    fragsRightReads = NULL
    for (i in 1:readLength) {
      fragsLeftPos = c(fragsLeftPos, fragsLeft[,1]+i-1)
      fragsLeftReads = c(fragsLeftReads, fragsLeft[,2])
      fragsRightPos = c(fragsRightPos, fragsRight[,1]+i-1)
      fragsRightReads = c(fragsRightReads, fragsRight[,2])
    }
    fragsLeft = cbind(fragsLeftPos, fragsLeftReads)
    fragsRight = cbind(fragsRightPos, fragsRightReads)
  }
  
  fragEnds = rbind(fragsLeft, fragsRight)
  fragEnds = fragEnds[order(fragEnds[,1]),]
  
  wigData = c(wigHeader, paste(fragEnds[,1], fragEnds[,2], sep = " "))
  
  return(wigData)
}


createVirtualFragmentLibraryMain <- function(totalFragments, totalFragmentsRev, firstCutter, secondCutter, readLength, onlyNonBlind = TRUE, useOnlyIndex = FALSE, minSize = 0, maxSize = -1, minFragEndSize = 0, maxFragEndSize = 10000000, chromosomeName = "chr1", libraryName = "default") {
  
  # if necessary, change chromosome name from "chr1"..."chrY" to "1"..."Y"
  if (useOnlyIndex) {
    totalFragments$chromosomeName = sub("chr", "", totalFragments$chromosomeName)
    totalFragmentsRev$chromosomeName = sub("chr", "", totalFragmentsRev$chromosomeName)
  }
  
  totalFragments["fragmentLength"] = nchar(totalFragments$fragmentSequences)
  totalFragmentsRev["fragmentLength"] = nchar(totalFragmentsRev$fragmentSequences)
  
  # calculate frag end lengths
  leftFragEndLength = ifelse(totalFragments$secondCutterFirst == -1, totalFragments$fragmentLength, totalFragments$secondCutterFirst - 1)
  rightFragEndLength = ifelse(totalFragments$secondCutterFirst == -1, totalFragments$fragmentLength, totalFragments$fragmentLength - totalFragments$secondCutterLast + 1 - nchar(secondCutter))
  leftFragEndLengthRev = ifelse(totalFragmentsRev$secondCutterFirst == -1, totalFragmentsRev$fragmentLength, totalFragmentsRev$secondCutterFirst - 1)
  rightFragEndLengthRev = ifelse(totalFragmentsRev$secondCutterFirst == -1, totalFragmentsRev$fragmentLength, totalFragmentsRev$fragmentLength - totalFragmentsRev$secondCutterLast + 1 - nchar(secondCutter))
  
  # valid 4C-seq reads map to one of the experiment's fragment ends and continue with the downstream sequence
  fragSeqStart = substr(totalFragments$fragmentSequences, start = 1, stop = readLength)
  fragSeqStartRev = substr(totalFragmentsRev$fragmentSequences, start = 1, stop = readLength)
  
  # check fragment ends for uniqueness
  uniqueFromStart = !duplicated(c(fragSeqStart, fragSeqStartRev))
  uniqueFromEnd = !duplicated(c(fragSeqStart, fragSeqStartRev), fromLast = TRUE)
  uniqueSeq = uniqueFromStart & uniqueFromEnd
  
  uniqueFragStart = uniqueSeq[1:nrow(totalFragments)]
  uniqueFragEnd = rev(uniqueSeq[(nrow(totalFragments)+1):(nrow(totalFragments)*2)])
  
  rm(totalFragmentsRev)
  
  # remove raw sequences
  totalFragments = totalFragments[,-5]
  
  # mark fragments with length < readLength as not usable, check minimum and maximum fragment end sizes
  minFragEndSize = max(readLength, minFragEndSize)
  leftLongEnough = ifelse((leftFragEndLength >= minFragEndSize & leftFragEndLength <= maxFragEndSize), TRUE, FALSE)
  rightLongEnough = ifelse((rightFragEndLength >= minFragEndSize & rightFragEndLength <= maxFragEndSize), TRUE, FALSE)
  
  leftFragEndValid = uniqueFragStart & leftLongEnough
  rightFragEndValid = uniqueFragEnd & rightLongEnough
  
  # centre of fragment: middle between first and last second cutter site, if second cutter present
  fragmentCentre = ifelse(totalFragments$secondCutterFirst == -1, round((totalFragments$fragmentStart + totalFragments$fragmentEnd) / 2), round((totalFragments$secondCutterFirst + totalFragments$secondCutterLast) / 2) + totalFragments$fragmentStart + 1)
  
  isNonBlind = totalFragments$secondCutterPresent
  
  fragData = data.frame(totalFragments[,1:3], fragmentCentre, isNonBlind, "fragmentLength" = totalFragments$fragmentLength, leftFragEndLength, rightFragEndLength, leftFragEndValid, rightFragEndValid)
  
  if (onlyNonBlind) {
    fragData = subset(fragData, fragData$isNonBlind == TRUE)
  }
  
  # if chosen, keep only fragments with a fragment length of more than X bp...
  fragData = subset(fragData, fragData$fragmentLength >= minSize)
  # ... and less than Y bp
  if (maxSize != -1) {
    fragData = subset(fragData, fragData$fragmentLength <= maxSize)
  }
  
  # output: fragment coordinates and flags for uniqueness
  if (libraryName == "default") {
    if (onlyNonBlind) {
      libraryName = paste("fragments_", firstCutter, "_", secondCutter, ".csv", sep = "")
    } else {
      libraryName = paste("fragments_", firstCutter, "_", secondCutter, "_with_blind_fragments", ".csv", sep = "")
    }
  }
  
  if (libraryName == "") {
    return(fragData)
  } else {
    write.table(fragData, file = libraryName, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  }
}



createVirtualFragmentLibrary_2 <- function(chosenGenome, firstCutter, secondCutter, readLength, onlyNonBlind = TRUE, useOnlyIndex = FALSE, minSize = 0, maxSize = -1, minFragEndSize = 0, maxFragEndSize = 10000000, useAllData = TRUE, chromosomeName = "chr1", libraryName = "default") {
  
  chromosomeNames = seqnames(chosenGenome)[grepl('^chr[0-9XY]{1,2}$', seqnames(chosenGenome))]
  
  totalFragments = NULL
  totalFragmentsRev = NULL
  
  # for each chromosome: split current chromosome at the first cutter sequence
  # and check for presence of the second cutter within the resulting fragments
  # --> remove non-unique and blind fragments (if chosen; default == TRUE) for final output
  for (i in 1:length(chromosomeNames)) {
    
    chromosomeToSplit = chosenGenome[[i]]
    
    if (class(chromosomeToSplit) == "MaskedDNAString") {
      chromosomeToSplit = unmasked(chromosomeToSplit)
    }
    
    if (useAllData) {
      message(paste("analyzing ", chromosomeNames[i], "...", sep = ""))
      chromosomeToSplitRev = reverseComplement(chromosomeToSplit)
      currentChromosome = splitChromosome(firstCutter, secondCutter, chromosomeToSplit, chromosomeNames[i])
      currentChromosomeRev = splitChromosome(firstCutter, secondCutter, chromosomeToSplitRev, chromosomeNames[i])
      totalFragments = rbind(totalFragments, currentChromosome)
      totalFragmentsRev = rbind(currentChromosomeRev, totalFragmentsRev)
    } else {
      
      if (length(chromosomeToSplit) > 20000000) {
        message(paste("analyzing ", chromosomeNames[i], "...", sep = ""))
        chromosomeToSplitRev = reverseComplement(chromosomeToSplit)
        currentChromosome = splitChromosome(firstCutter, secondCutter, chromosomeToSplit, chromosomeNames[i])
        currentChromosomeRev = splitChromosome(firstCutter, secondCutter, chromosomeToSplitRev, chromosomeNames[i])
        totalFragments = rbind(totalFragments, currentChromosome)
        totalFragmentsRev = rbind(currentChromosomeRev, totalFragmentsRev)
      } else {
        message(paste("skipping ", chromosomeNames[i], " due to its length", sep = ""))
      }
    }
  }
  
  if (!(is.null(totalFragments))) {
    createVirtualFragmentLibraryMain(totalFragments, totalFragmentsRev, firstCutter, secondCutter, readLength, onlyNonBlind, useOnlyIndex, minSize, maxSize, minFragEndSize, maxFragEndSize, chromosomeName, libraryName)
  } else {
    message("No fragments created; please use 'useAllData = TRUE' for genomes with small chromosomes")
  }
}


splitChromosome <- function(firstCutter, secondCutter, chromosomeToSplit, chromosomeName) {
  
  # first step: get position of the cutter sequences and calculate start and end of fragment in between
  #  --> cutter sequence neglected for fragment to provide disjunct fragment intervals
  rawFragments = matchPattern(firstCutter, chromosomeToSplit)
  fragmentStart = c(1, end(rawFragments) + 1)
  fragmentEnd = c(start(rawFragments) - 1, length(chromosomeToSplit))
  
  # second step: read chromosome as string and check which fragments are blind (no second cutter present
  # or non-blind (second cutter sequence present within fragment sequence)
  chromosomeTotal = toString(chromosomeToSplit)
  fragmentSequences = unlist(strsplit(chromosomeTotal, split=toupper(firstCutter)))
  secondCutterPresent = grepl(toupper(secondCutter), fragmentSequences)
  
  # sequences at start and end of chromosome not counted as non-blind fragments 
  # --> valid non-blind fragments: only [FCE ... SCE ... FCE], not [1 ... (SCE) ... FCE] or [FCE ... (SCE) ... END] 
  secondCutterPresent[1] = FALSE
  secondCutterPresent[length(secondCutterPresent)] = FALSE
  
  emptyLastFrag = FALSE
  
  if (length(fragmentEnd) > length(secondCutterPresent)) {
    secondCutterPresent = c(secondCutterPresent, FALSE)
    fragmentSequences = c(fragmentSequences, "")
    emptyLastFrag = TRUE
  }
  
  # fragments total
  fragmentTable = data.frame(chromosomeName, fragmentStart, fragmentEnd, secondCutterPresent, fragmentSequences)
  fragmentTable[,5] = as.vector(fragmentTable[,5])
  
  # delete possible empty fragment if cutter sequence is at the end of the genome
  if (emptyLastFrag) {
    fragmentTable = fragmentTable[-nrow(fragmentTable),]
  }
  
  secondCutterPos = gregexpr(toupper(secondCutter), fragmentTable[,5])
  secondCutterFirst = NULL
  secondCutterLast = NULL
  
  for (i in 1:length(secondCutterPos)) {
    secondCutterFirst[i] = secondCutterPos[[i]][1]
    secondCutterLast[i] = secondCutterPos[[i]][length(secondCutterPos[[i]])]
  }
  
  fragmentTable["secondCutterFirst"] = secondCutterFirst
  fragmentTable["secondCutterLast"] = secondCutterLast
  
  # return bed-like format: chromosome name, fragment start, fragment end plus fragment sequence and second cutter positions
  return(fragmentTable)
}


readsToFragments_2 <- function(expData, fragmentLib) {
  
  readsTotal = rawReads(expData)
  # check all chromosomes where reads are present
  chromosomeNames = unique(as.vector(seqnames(readsTotal)))  # why is chr missing here?
  chromosomeNames = chromosomeNames[grepl('^(chr([1-9]|1[0-9]|2[0-1]|X|Y)|([1-9]|1[0-9]|2[0-1]|X|Y))$', chromosomeNames)]  
  #if(!unique(grepl('chr', chromosomeNames))){  # TODO check this and the upper regex, rm when checked chr before
  #  chromosomeNames <- paste0('chr', chromosomeNames)
  #}
  # read fragment data (with  meta data, like uniqueness of frag-ends + second cutter site)
  fragmentTableTotal = read.table(fragmentLib, header = TRUE)
  
  fragmentDataTotal = NULL
  
  for (i in 1:length(chromosomeNames)) {
    # message(chromosomeNames[i])
    reads = subset(readsTotal, as.character(seqnames(readsTotal)) == chromosomeNames[i])
    
    readsPlus = subset(reads, as.character(strand(reads)) == "+")
    readsMinus = subset(reads, as.character(strand(reads)) == "-")
    
    if (length(readsPlus) > 0 & length(readsMinus) > 0) {
      
      baseCovPlus = as.numeric(coverage(readsPlus)[[chromosomeNames[i]]])
      baseCovMinus = as.numeric(coverage(readsMinus)[[chromosomeNames[i]]])
      
      if (startsWith(x = chromosomeNames[i], 'chr')){   # TODO rm if corrected chr earlier
        fragmentTable = subset(fragmentTableTotal, fragmentTableTotal$chromosomeName == chromosomeNames[i])
      } else {
       fragmentTable = subset(fragmentTableTotal, fragmentTableTotal$chromosomeName == paste0('chr', chromosomeNames[i]))
      }
      
      # add zeroes if last fragment is not covered
      tempPlus = fragmentTable[nrow(fragmentTable),3] - length(baseCovPlus)
      if (tempPlus < 0) {
        tempPlus = 0
      }
      baseCovPlus = c(baseCovPlus,  rep(0, times=tempPlus))
      tempMinus = fragmentTable[nrow(fragmentTable),3] - length(baseCovMinus)
      if (tempMinus < 0) {
        tempMinus = 0
      }
      baseCovMinus = c(baseCovMinus,  rep(0, times=tempMinus))
      
      if (fragmentTable$fragmentEnd[1] == 0) {
        fragmentTable$fragmentEnd[1] = 1
      }
      
      # calculate read count for fragment start, end, and average
      leftFragEndReads = baseCovPlus[fragmentTable$fragmentStart]
      rightFragEndReads = baseCovMinus[fragmentTable$fragmentEnd]
      fragEndReadsAverage = round((leftFragEndReads + rightFragEndReads) / 2)
      
      fragmentData = data.frame(fragmentTable, leftFragEndReads, rightFragEndReads, fragEndReadsAverage)
      fragmentDataTotal = rbind(fragmentDataTotal, fragmentData)
    }
  }
  
  return(fragmentDataTotal)  
}

# # Vis VP ####
# 
# visualizeViewpointMain <- function(fragsToVisualize, poi = data.frame(chr = character(), start = character(), end = character(), name = character(), colour = character()), plotFileName = "", windowLength = 5, interpolationType = "median", picDim = c(9, 5), maxY = -1, minQuantile = 0.2, maxQuantile = 0.8, mainColour = "blue", plotTitle = "4C-seq plot", loessSpan = 0.1, xAxisIntervalLength = 50000, yAxisIntervalLength = 500, useFragEnds = TRUE) {
# 
#   # prepare different color shades
#   colours = c(paste(mainColour, "", sep = ""), paste(mainColour, "2", sep = ""), paste(mainColour, "3", sep = ""))
# 
#   # prepare data
#   yVal = 0
# 
#   position = round((fragsToVisualize$start + fragsToVisualize$end) / 2)
#   averageReads = fragsToVisualize$reads
# 
#   if (interpolationType == "median") {
#     # smoothed data (median with defined window length)
#     yVal = runmed(averageReads, windowLength, endrule = "constant")
# 
#   } else if (interpolationType == "mean") {
#     # smoothed data (mean with defined window length)
#     yVal = runmean(averageReads, windowLength, endrule = "constant")
# 
#   } else {
#     # raw data
#     yVal = averageReads
#   }
# 
#   # if no (reasonable) maximum y-value is provided, use the maximum of yValA and yValB, rounded up to the next y-axis marker
#   if (maxY == -1) {
#     maxY = ceiling(max(yVal) / 500) * 500
#   }
# 
#   # prepare plot
#   if (plotFileName != "") {
#     if (grepl(".pdf", plotFileName)) {
#       pdf(file = plotFileName, title = plotTitle, width = picDim[1], height = picDim[2], useDingbats = FALSE)
#     } else {
#       tiff(filename = plotFileName, width = picDim[1], height = picDim[2], compression = "none", bg = "white", pointsize = 20)
#     }
#   }
# 
#   # basic plot
#   plot(position, yVal,  ylim = c(0,maxY), ylab = "fragment read count (RPM)", xlab = "fragment position", type = "n", main = plotTitle, lty = 1, axes = FALSE)
# 
#   # plot sample data
#   basic4CPlot(position, yVal, averageReads, minQuantile, maxQuantile, windowLength, colours, loessSpan)
# 
#   drawMetaData(fragsToVisualize$start[1], fragsToVisualize$end[nrow(fragsToVisualize)], maxY, poi, xAxisIntervalLength, yAxisIntervalLength)
# 
#   if (plotFileName != "") {
#     dev.off()
#   }
# 
#   # prepare legend plot
#   scale = 100
# 
#   if (plotFileName != "") {
#     if (grepl(".pdf", plotFileName)) {
#       pdf(file = "quantile_legend.pdf", title = "Main Trend", width = 2, height = 5, useDingbats = FALSE)
#     } else {
#       tiff(filename = "quantile_legend.tiff", width = 200, height = 500, compression = "none", bg = "white", pointsize = 20)
#     }
# 
#     # basic plot
#     plot(c(0,10), c(0,1), type='n', main="Main trend", axes = FALSE, xlab = '', ylab = "loess-smoothed quantiles")
# 
#     # add axis...
#     axis(2, c(0, 0.5, 1), labels = c(paste(minQuantile, "%", sep = ""), "50%", paste(maxQuantile, "%", sep = "")), las = 1)
# 
#     # ... quantile representation...
#     for (i in 1:100) {
#       y = (i-1) / scale
#       rect(0, 0, 10, 1, col = "grey90", border=NA)
#     }
# 
#     # ... and trend line
#     lines(c(0, 10), c(0.5, 0.5), col = "black", lwd = 3)
# 
#     dev.off()
#   }
# }
# visualizeViewpoint <- function(expData, poi = data.frame(chr = character(), start = character(), end = character(), name = character(), colour = character()), plotFileName = "", windowLength = 5, interpolationType = "median", picDim = c(9, 5), maxY = -1, minQuantile = 0.2, maxQuantile = 0.8, mainColour = "blue", plotTitle = "4C-seq plot", loessSpan = 0.1, xAxisIntervalLength = 50000, yAxisIntervalLength = 500, useFragEnds = TRUE) {
# 
#   fragsToVisualize = formatFragmentData(expData = expData, useFragEnds)
#   poi = pointsOfInterest(expData)
# 
#   visualizeViewpointMain(fragsToVisualize, poi, plotFileName, windowLength, interpolationType, picDim, maxY, minQuantile, maxQuantile, mainColour, plotTitle, loessSpan, xAxisIntervalLength, yAxisIntervalLength, useFragEnds)
# }
# 
# drawMetaData <- function(minX, maxX, maxY, poi = data.frame(chr = character(), start = character(), end = character(), name = character(), colour = character()), xAxisIntervalLength = 50000, yAxisIntervalLength = 500) {
# 
#   # if present, print points of interest as defined by the user
#   if (nrow(poi) > 0) {
#     pointsX = c(poi$start, poi$end)
#     pointsY = round(maxY * 0.97)
# 
#     points(round((poi$start + poi$end) / 2), rep(round(maxY * 0.95), times = nrow(poi)), cex = 1.0, pch = 25, col = as.vector(poi$colour), bg = as.vector(poi$colour))
#     text(round((poi$start + poi$end) / 2), round(maxY * 0.99), poi$name, cex = 1.0)
#   }
# 
#   # print axes with custom intervals
#   minXAxis = floor(minX / xAxisIntervalLength) * xAxisIntervalLength
#   maxXAxis = ceiling(maxX / xAxisIntervalLength) * xAxisIntervalLength
#   seqXAxis = seq(minXAxis, maxXAxis, by = xAxisIntervalLength)
#   axis(side = 1, at = seqXAxis, labels = seqXAxis, las = 0)
# 
#   seqYAxis = seq(0, ceiling(maxY / yAxisIntervalLength) * yAxisIntervalLength, by = yAxisIntervalLength)
#   axis(side = 2, at = seqYAxis, labels = seqYAxis, las = 2)
# }
# 

# # Format Fragment Data ####
# formatFragmentData <- function(expData, useFragEnds = TRUE) {
# 
#   normalizedFragData = nearCisFragments(expData)
# 
#   formatFragmentDataMain(normalizedFragData, useFragEnds)
# }


# .formatFragmentData_DataFrame <- function(expData, useFragEnds = TRUE) {
#
#   normalizedFragData = expData
#
#   formatFragmentDataMain(normalizedFragData, useFragEnds)
# }
# 
# 
# formatFragmentDataMain <- function(normalizedFragData, useFragEnds = TRUE) {
# 
#   if (useFragEnds) {
#     fragmentDataLeftUnique = subset(normalizedFragData, normalizedFragData$leftFragEndValid == TRUE)
#     fragmentDataRightUnique = subset(normalizedFragData, normalizedFragData$rightFragEndValid == TRUE)
# 
#     fragsLeft = fragmentDataLeftUnique[,c(1,2,4,11)]
#     fragsRight = fragmentDataRightUnique[,c(1,4,3,12)]
#     fragsRight[,2] = fragsRight[,2] + 1
# 
#     colnames(fragsLeft) = c("chrom", "start", "end", "reads")
#     colnames(fragsRight) = c("chrom", "start", "end", "reads")
# 
#     fragEnds = rbind(fragsLeft, fragsRight)
#     fragmentDataFinal = fragEnds[order(fragEnds[,1], fragEnds[,2]),]
# 
#   } else {
#     fragments = subset(normalizedFragData, normalizedFragData$leftFragEndValid == TRUE & normalizedFragData$rightFragEndValid == TRUE)
# 
#     fragmentDataFinal = fragments[,c(1:3,13)]
#     colnames(fragmentDataFinal) = c("chrom", "start", "end", "reads")
#   }
# 
#   row.names(fragmentDataFinal) = NULL
# 
#   return(fragmentDataFinal)
# }
# 
# 
# 
# 
# 
# # DrawHeatmap ####
# drawHeatmap_Data4Cseq <- function(expData, plotFileName = "", smoothingType = "median", picDim = c(9, 2.2), bands = 5, cutoffLog = -7.0, xAxisIntervalLength = 50000, legendLabels = expression(2^-7, 2^0), useFragEnds = TRUE) {
# 
#   fragsToVisualize = formatFragmentData(expData, useFragEnds)
# 
#   drawHeatmapMain(fragsToVisualize, plotFileName, smoothingType, picDim, bands, cutoffLog, xAxisIntervalLength, legendLabels, useFragEnds)
# }
# 
# 
# drawHeatmap_DataFrame <- function(expData, plotFileName = "", smoothingType = "median", picDim = c(9, 2.2), bands = 5, cutoffLog = -7.0, xAxisIntervalLength = 50000, legendLabels = expression(2^-7, 2^0), useFragEnds = TRUE) {
#   
#   fragsToVisualize = expData
#   
#   if (!useFragEnds) {
#     print("drawHeatmap: use of fragment-wise interpolated data not possible for data frame input")
#   }
#   
#   drawHeatmapMain(fragsToVisualize, plotFileName, smoothingType, picDim, bands, cutoffLog, xAxisIntervalLength, legendLabels, useFragEnds)
# }
# 
# 
# drawHeatmapMain <- function(fragsToVisualize, plotFileName = "", smoothingType = "median", picDim = c(9, 2.2), bands = 5, cutoffLog = -7.0, xAxisIntervalLength = 50000, legendLabels = expression(2^-7, 2^0), useFragEnds = TRUE) {
#   
#   position = round((fragsToVisualize$start + fragsToVisualize$end) / 2)
#   averageReads = fragsToVisualize$reads
#   
#   # prepare heatmap-like data (varied window lengths for running mean / running median)
#   signal = matrix(0, bands, nrow(fragsToVisualize))
#   
#   for (i in 1:bands) {
#     
#     wl = i * 2 - 1
#     
#     if (smoothingType == "median") {
#       # smoothed data (median with defined window length)
#       signal[i,] = runmed(averageReads, wl, endrule = "keep")
#       
#     } else {
#       # smoothed data (mean with defined window length)
#       signal[i,] = runmean(averageReads, wl, endrule = "keep")
#     }     
#   }
#   
#   # prepare heatmap plot
#   if (plotFileName != "") {
#     if (grepl(".pdf", plotFileName)) {
#       pdf(file = plotFileName, title = plotFileName, width = picDim[1], height = picDim[2], useDingbats = FALSE)  
#     } else {
#       tiff(filename = plotFileName, width = picDim[1], height = picDim[2], compression = "none", bg = "white", pointsize = 20)
#     }
#   }
#   
#   plot(position, averageReads,  ylim = c(0,1), ylab = "wl", xlab = "fragment position", type = "n", lty = 1, axes = FALSE)
#   
#   scale = 1 / (bands + 1.5)
#   maxSignal = max(signal)
#   
#   # log-scale signal
#   signal = log2(signal/maxSignal)
#   
#   # set cut-offs for log-scaled data
#   capMin = cutoffLog
#   capMax = 0.0
#   
#   # define axes with custom intervals
#   minXAxis = floor(fragsToVisualize$start[1] / xAxisIntervalLength) * xAxisIntervalLength
#   maxXAxis = ceiling(fragsToVisualize$end[nrow(fragsToVisualize)] / xAxisIntervalLength) * xAxisIntervalLength
#   seqXAxis = seq(minXAxis, maxXAxis, by = xAxisIntervalLength)
#   axis(side = 1, at = seqXAxis, labels = seqXAxis, las = 0)
#   
#   seqYAxis = c(0.5, bands - 0.5) * scale
#   labelsYAxis = c(2*bands-1,1)    
#   
#   axis(side = 2, at = seqYAxis, labels = labelsYAxis, las = 2)
#   
#   # mark unused fragments (blind / repeated)
#   rect(minXAxis, 0, maxXAxis, bands * scale, col = "black", border = NA)
#   
#   # prepare color palette 
#   colorNumber = capMax*100-capMin*100 + 1
#   colors = heat.colors(colorNumber)
#   
#   # pick colours for the fragments
#   for (i in 1:bands) {
#     
#     for (j in 1:nrow(fragsToVisualize)) {
#       
#       if (signal[i,j] > capMax) {
#         signal[i,j] = capMax
#       }
#       if (signal[i,j] < capMin) {
#         signal[i,j] = capMin
#       }
#       
#       index = round((signal[i,j] - capMin) * 99) + 1
#       
#       pickedColor = colors[index]
#       
#       # mark 'missing' fragments at the sides (only fully present running median / mean windows are used)
#       if ((j < i) | (j > (nrow(fragsToVisualize) - i))) {
#         pickedColor = "black"
#       }
#       
#       # print fragments in chosen colours
#       rect(fragsToVisualize[,2][j], (bands-i+1)*scale, fragsToVisualize[,3][j], (bands-i)*scale, col = pickedColor, border = NA)
#       
#     }
#   }
#   
#   if (plotFileName != "") {
#     dev.off()
#   }
#   
#   # prepare heatmap legend
#   chosenPalette = heat.colors(100)
#   
#   # prepare color legend plot
#   if (plotFileName != "") {
#     if (grepl(".pdf", plotFileName)) {
#       pdf(file = "color_legend.pdf", title = "Color Legend", width = 2, height = 5, useDingbats=FALSE)  
#     } else {
#       tiff(filename = "color_legend.tiff", width = 200, height = 500, compression = "none", bg = "white", pointsize = 20)
#     }
#     
#     par(mar=c(5.1,4.1,4.1,4.1))
#     
#     # basic plot
#     plot(c(0,10), c(0,1), type='n', main="Signal intensity", axes = FALSE, xlab = "", ylab = "relative window coverage (log-scaled)")
#     
#     # add axis...
#     legendLength = length(legendLabels) - 1
#     axis(2, c((0:legendLength)/legendLength), labels = legendLabels, las = 1)
#     
#     # ... and color bands
#     for (i in 1:(length(chosenPalette)-1)) {    
#       y = (i-1) / 100
#       rect(0, y, 10, y+1/100, col = chosenPalette[i], border=NA)
#     }
#     
#     dev.off()
#   }
# }
# 
# 
# 
# 
# setMethod("drawHeatmap",
#           signature=signature(expData="Data4Cseq"),
#           drawHeatmap_Data4Cseq)
# 
# setMethod("drawHeatmap",
#           signature=signature(expData="data.frame"),
#           drawHeatmap_DataFrame)
# 
# setMethod("drawHeatmapMain",
#           signature=signature(fragsToVisualize="data.frame"),
#           drawHeatmapMain)

