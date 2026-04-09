#########################################
##
## fourSig
## Author: Joshua Starmer <josh.starmer@gmail.com>, 2013
## 
## Copyright (C) 2013, Joshua Starmer
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License
## as published by the Free Software Foundation; either version 2
## of the License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
##
## NOTE: This code was based on a program written
## by Elzo de Wit and I am very much in his debt to him for sending me his
## work.  The modifications and additions that I have made are so substantial
## that this fork deserves a new name, so as not to be confused with his
## program.  However, I have retained one of his "helper" functions because
## it is useful and efficient.  Any code that retains his original
## algorithm has been annotated as souch.  Also, any clever R tricks in this
## code were inspired by Elzo's work.
##
#########################################

source("readsToWindows.R")
source("plotReads.R")

#########################################################
##
## GLOBAL CONSTANTS THAT MAKE EVERYTHING WAY EASIER TO READ
##
#########################################################


TAB.CHR        <- 1
TAB.POS        <- 2
TAB.READ.COUNT <- 3
TAB.MAP        <- 4
TAB.4BP        <- 5
TAB.SHORT      <- 6
TAB.6BP        <- 7
TAB.LENGTH     <- 8

DATA.CHR        <- 1
DATA.POS        <- 2
DATA.READ.COUNT <- 3
DATA.MAP        <- 4
DATA.4BP        <- 5
DATA.SHORT      <- 6
DATA.6BP        <- 7
DATA.LENGTH     <- 8

NUM.ITERATIONS   <- 100
FDR.CUTOFF       <- 0.01
FDR.PROB         <- 0.05



#########################################################
##
## FUNCTIONS
##
#######################################################

## runningSum, by Elzo de Wit,
## Quickly calculate the running sum
## Splinter et al. 2012 Methods
## Supplemental R script
runningSum <- function( x, n ) {
  sum.v <- c(0,cumsum(x))

  ## diff, in this case is a vector of
  ## sum.v[(1+n):length(sum.v)] - sum.v[1:(length(sum.v)-n)]
  ## i.e.
  ## sum.v = 0 0 1 1 1 2 2 2 3 3 3
  ##         ^       ^               1 - 0 = 1
  ## sum.v = 0 0 1 1 1 2 2 2 3 3 3
  ##           ^       ^             2 - 0 = 2
  ## sum.v = 0 0 1 1 1 2 2 2 3 3 3
  ##             ^       ^           2 - 1 = 1 etc...
  diff(sum.v,n)
}	


shuffleRaw <- function(read.data) {
  total.reads = sum(read.data)
  num.bins = length(read.data)
  
  new.data = round(runif(min=1, max=num.bins, n=total.reads))
  return(table(factor(new.data, levels=1:num.bins))[])
}


fdrShuffleMemoryless <- function(read.data, window=20,
                                 iterations=NUM.ITERATIONS,
                                 fdr, fdr.prob,
                                 raw=FALSE) {
  
  max.level = sum(head(sort(read.data, decreasing=TRUE), n=window))
  
  real.data.table <- table(factor(runningSum(read.data, window),
                                  levels=0:max.level) )
  
  max.reads.in.bin <-
    max(as.numeric(names(real.data.table[real.data.table > 0])))
  
  real.data.sum <- cumsum(rev(real.data.table))
  
  ## now generate a table of randomized data.
  random.data.table  <- NULL
  shuffled.min.reads <- NULL
  
  for(i in 1:iterations) {
    random.data.table <- NULL

    random.reads = shuffleRaw(read.data)
    random.data.table <-
      table(factor(runningSum(random.reads,window),levels=0:max.level))
    
    random.data.sum <- cumsum(rev(random.data.table))

    fractions <- random.data.sum[real.data.sum > 0] /
      real.data.sum[real.data.sum > 0]

    sig <- fractions[fractions < fdr]
  
    min.read <- NaN
    if (length(sig) > 0) {
      min.read <- min(as.numeric(names(sig)))
    }

    shuffled.min.reads[i] <- min.read
  }

  possible.levels <- levels(as.factor(shuffled.min.reads))

  print(paste("  These values were proposed for FDR cut offs ( max",
              max.reads.in.bin, "):", possible.levels))

  level.counts <- c()
  fdr.to.use   <- tail(possible.levels, n=1) # start with the largest value
  prob.under.estimate <- 1
  largest.fdr.found   <- FALSE;
  for(i in possible.levels) {
    level.count   <- length(shuffled.min.reads[shuffled.min.reads == i])
    percent.level <- (level.count / iterations)

    print(paste(" ", i, "occurred in", round((percent.level*100), digits=2),
                "percent of", iterations, "shufffles"))
    
    print(paste("  prob.under.estimate:", round(prob.under.estimate, digits=2)))
    
    if ((!largest.fdr.found) && (prob.under.estimate <= fdr.prob)) {
      largest.fdr.found <- TRUE
      print(paste("    setting fdr.to.use to", i))
      fdr.to.use <- i
    }

    prob.under.estimate <- prob.under.estimate - percent.level

    level.counts <- c(level.counts, i, level.count)

  }
  if (!largest.fdr.found) {
    fdr.to.use <- as.numeric(fdr.to.use) + 1
  }

  print(paste("  Using", fdr.to.use, "for the cutoff"))
  
  return (list(max.reads.real=max.reads.in.bin,
               fdr.cut.off=fdr.to.use,
               fdr.hist=level.counts))
}

################################
##
## BEGIN: getDomains()
##
################################
getDomains <- function(full.data, shuffle.data, 
                       window=500,
                       iterations=NUM.ITERATIONS,
                       raw=FALSE,
                       fdr, fdr.prob,
                       fixed.cutoff=NULL, min.chr.read.count=0) {  

  domains     <- NULL
  chromosomes <- unique(full.data[,DATA.CHR])
  
  print("")
  
  for( chr in chromosomes ) {
    chr.data = data.frame(chr=full.data[full.data[,DATA.CHR]==chr,DATA.CHR],
      pos       = full.data[full.data[,DATA.CHR]==chr,DATA.POS],
      read.count= full.data[full.data[,DATA.CHR]==chr,DATA.READ.COUNT],
      map       = full.data[full.data[,DATA.CHR]==chr,DATA.MAP])
      ##cut.4bp   = full.data[full.data[,DATA.CHR]==chr,DATA.4BP],
      ##short     = full.data[full.data[,DATA.CHR]==chr,DATA.SHORT],
      ##cut.6bp   = full.data[full.data[,DATA.CHR]==chr,DATA.6BP],
      ##length    = full.data[full.data[,DATA.CHR]==chr,DATA.LENGTH])
      
    chr.data.shuffle = shuffle.data[shuffle.data[,DATA.CHR]==chr,]

    chr.read.count <- sum(chr.data[,DATA.READ.COUNT])
    
    print(paste("Working on chr", chr, sep=""))
    print(paste("Total reads on chr", chr , ": ", chr.read.count, sep=""))
    
    if (chr.read.count < min.chr.read.count) {
      print(paste("Skipping chr", chr, " since it has fewer than ", min.chr.read.count, " reads mapped to it", sep=""))
      print("");
      next;
    }
    
    ## make sure there are more fragments than the window size.
    num.frags <- nrow(chr.data)
    if (nrow(chr.data) < window) {
      print(paste("Warning: the window size, ", window, ", is larger than the number of fragments, ", num.frags, sep=""))
      print(paste("Skipping chromosome ", chr))
      next;
    }
    
    ## if there are more fragments than the window size, but still not many
    ## fragments, let the user know, but don't quit
    if(num.frags < 1000) {
      print(paste("Warning:", chr, "has less than 1000 fragments"))
      print(paste(".... Attempting to generate stats for ",
                  num.frags,
                  "fragments..."))
    }
    
    print(paste("Shuffling", chr))

    read.data.shuffle <- chr.data.shuffle[,DATA.READ.COUNT]
    read.data.full    <- chr.data[,DATA.READ.COUNT]
    
    pos <- chr.data[, DATA.POS]

    min.reads.for.fdr <- NULL
    if (is.null(fixed.cutoff)) {
      min.reads.for.fdr <- fdrShuffleMemoryless(read.data.shuffle,
                                                window=window,
                                                raw=raw,
                                                iterations=iterations,
                                                fdr=fdr,
                                                fdr.prob=fdr.prob)
      
      cut.off <- as.numeric(min.reads.for.fdr$fdr.cut.off)
      
      if(is.nan(as.numeric(cut.off))) {
        print(paste("Chromosome", chr, "had no enriched regions"))
        print("")
        next
      } else {
        print("")
      }
      
    } else {      
      cut.off <- as.numeric(fixed.cutoff)
    }
    
    
    ## pick out the positions that have more reads in them
    window.data <- runningSum(read.data.full, window)
    sig.indices <- which(window.data >= cut.off)
    if (length(sig.indices) == 0) {
      print(paste("There were no windows that had more reads than", cut.off))
      print("")
    } else {

      ##print(paste("length(sig.indices): ", length(sig.indices), sep=""))
      ##print(paste("sig.indices[1]: ", sig.indices[1], sep=""))
      ##print(paste("read count: ", window.data[sig.indices[1]], sep=""))
      
      extra.sig <- matrix(nrow=length(sig.indices), ncol=2)

      ## here we test to see if the region/window is sigificant because
      ## of PCR amplification of a single fragment.  We do this by
      ## removing the fragment with the most reads (or replacing that
      ## fragment with the average of the neighboring fragments) and
      ## see if the window is still sigificant.
      peak.counter <- 1
      
      for (enriched.i in sig.indices) {
        fragment.reads <- read.data.full[enriched.i:(enriched.i + (window-1))]
        total.reads    <- sum(fragment.reads)
        
        if (is.na(total.reads)) {
          if (is.na(fragment.reads)) {
            print("fragment.reads also had no value")
          } else {
            print(paste("fragment.reads", fragment.reads))
          }
          
          print(paste("enriched.i", enriched.i))
          print(paste("window", window))
          new.dim = dim(chr.data)
          print(paste("new.dim", new.dim))
        }
        
        max.pos     <- 0
        max.count   <- 0
        left.count  <- 0
        right.count <- 0
        for(i in 1:window) {
          if (fragment.reads[i] > max.count) {
            max.count <- fragment.reads[i]
            max.pos   <- i
            if (i > 1) {
              left.count <- fragment.reads[(i-1)]
            }
            if (i < window) {
              right.count <- fragment.reads[(i+1)]
            }
          }
        }
        
        avg.count = (left.count + right.count) / 2
        
        if (is.na(total.reads)) {
          print("total reads had no value!")
        }
        if (is.na(max.count)) {
          print("max count had no values!")
        }
        if (is.na(cut.off)) {
          print("cut.off had no value!")
        }
        
        if ((total.reads - max.count) >= cut.off) {
          extra.sig[peak.counter,] <- c(enriched.i, 1)
        } else if ((total.reads - max.count + avg.count) >= cut.off) {
          extra.sig[peak.counter,] <- c(enriched.i, 2)
        } else {
          extra.sig[peak.counter,] <- c(enriched.i, 3)
        }
        peak.counter <- peak.counter + 1
      } # END for

      ##print(paste("peak.counter = ", peak.counter, sep=""))
      ## now convert the enriched positions into domains
      start.indices <- vector()
      end.indices   <- vector()
      
      ## first find the first and last window for each domain
      print("Determining peak coordinates")
      index.i <- 1
      start.end.index <- 1
      ##print(paste("sig.indices:", sig.indices))
      diff.vector <- diff(x=sig.indices, lag=1)
      ##print(paste("diff.vector:", diff.vector))
      diff.vector.length <- length(diff.vector)
      ##print(paste("diff.vector.length:", diff.vector.length))
      
      if (diff.vector.length == 0) {
        start.indices[start.end.index] <- sig.indices[1] 
        end.indices[start.end.index] <- sig.indices[1]
      } else {
        while(index.i <= diff.vector.length) {
          start.indices[start.end.index] <- sig.indices[index.i]
          
          done <- FALSE
          while(!done & (index.i <= diff.vector.length)) {
            if (diff.vector[index.i] == 1) {
              index.i <- index.i + 1
            } else {
              done <- TRUE
              end.indices[start.end.index] <- sig.indices[index.i]
            }
          }
          if (!done) {
            end.indices[start.end.index] <- sig.indices[index.i]
          }
          start.end.index <- start.end.index + 1
          index.i <- index.i + 1
        }
      }
      #print(paste("length(start.indices): ", length(start.indices), sep=""))
      #print(paste("start.indices: ", start.indices, sep=""))
      
      
      print("Trimming windows")
      start.end.index <- 1
      while(start.end.index <= length(start.indices)) {
        start.frag.index <- start.indices[start.end.index]

        ## make sure the current start index is pointing to a position that
        ## has reads.... (trim the 5' end)
        found.reads <- FALSE
        while(!found.reads & (start.frag.index <= num.frags)) {
          if (read.data.full[start.frag.index] > 0) {
            found.reads <- TRUE
          } else {
            start.frag.index <- start.frag.index + 1
          }
        }
        if (!found.reads) {
          cat("ERROR: There were no reads after a significant peak start!\n")
          return(NULL)
        }
        start.indices[start.end.index] <- start.frag.index

        ## now figure out where the peak ends...
        end.frag.index <- (end.indices[start.end.index] + window - 1)
        if (end.frag.index > num.frags) {
          end.frag.index <- num.frags
        }
        
        found.reads <- FALSE
        while(!found.reads & (end.frag.index >= start.frag.index)) {
          if (read.data.full[end.frag.index] > 0) {
            found.reads = TRUE
          } else {
            end.frag.index <- end.frag.index - 1
          }
        }
        end.indices[start.end.index] <- end.frag.index
        
        start.end.index <- start.end.index + 1

      } # END while(start.end.index <= length(start.indices)) {
      
      ## convert the 3C fragment indices into base pair coordinates
      print("Coverting fragment indices to base pair coordinates")
      pos <- chr.data[, DATA.POS]
      extra.sig[,1] <- pos[as.numeric(extra.sig[,1])]
      start <- pos[start.indices]
      end   <- pos[end.indices]
      
      chr.dom <- list(start=start, end=end)
      if (length(chr.dom$start) == 0) {
        print(paste("Chr", chr, " did not have any enriched domains", sep=""))
        print("");
        next;
      } else {
        #print(paste("Chr", chr, " had ", length(chr.dom$start), " domains", sep="")) 
      }
      
      sig.vals <- NULL
      
      max.sig.pos <- extra.sig[nrow(extra.sig),1]
      
      for(i in 1:length(chr.dom$start)) {
        
        if ((chr.dom$start[i] <= max.sig.pos) &
            (chr.dom$end[i] <= max.sig.pos)) {
          extra.sig.in.domain <- 
            extra.sig[((extra.sig[,1] >= chr.dom$start[i])
                       &
                       (extra.sig[,1] <= chr.dom$end[i])),]
          
          if (is.vector(extra.sig.in.domain)) {
            sig.vals[i] <- extra.sig.in.domain[2]
          } else {
            sig.vals[i] <- min(extra.sig.in.domain[,2])
          }
        } else {
          sig.vals[i] <- extra.sig[nrow(extra.sig), 2]
        }
      } # END for
        

      sig.counts <- NULL
      for(i in 1: length(chr.dom$start)) {              
        sig.counts[i] <- sum(chr.data[chr.data$pos >= chr.dom$start[i] &
                                      chr.data$pos <= chr.dom$end[i],]$read.count)
      }
      
      ## package everything up
      chr.dom <- data.frame(chr   = rep(chr, length(chr.dom$start)),
                            start = chr.dom$start,
                            end   = chr.dom$end,
                            category= sig.vals,
                            reads = sig.counts)
      
      ## add this chromosome's domains to the other domains
      domains <- rbind(domains, chr.dom, deparse.level=0)
    } # END } else { - the else is what you do if there are significant regions
  }
  domains
}	
################################
##
## END: getDomains()
##
################################



readTab <- function( fileName, onlyMappable=TRUE){
  tab.out <- read.table(file=fileName, sep="\t", header=TRUE)
  
  ## now remove all of the fragments that can't be mapped...
  mappable <- NULL
  if (onlyMappable) {
    print("Only using fragments that have the mappable flag")
    mappable <- tab.out[tab.out[TAB.MAP] == TRUE,]
  } else {
    mappable <- tab.out
  }
  mappable
}


## vp = viewpoint - the known part of 4c
## this nomenclature is from Splinter et al. 2012
fourSig <- function(filename, chr,
                    auto.mask = NULL,
                    mask.start=NULL,
                    mask.end=NULL,
                    invert.mask=FALSE,
                    out.file=NULL,
                    cis.method="raw",
                    trans.method="raw",
                    trans.only=FALSE,
                    cis.only=FALSE,
                    window.size=20,
                    iterations=NUM.ITERATIONS,
                    fdr=FDR.CUTOFF,
                    fdr.prob=FDR.PROB,
                    only.mappable=TRUE,
                    fixed.cutoff=NULL,
                    print.results=TRUE) { 

  print(paste("Reading in", filename))
  data <- readTab(filename, onlyMappable=only.mappable) # read in the data...
  print("Done reading the file");

  if (!trans.only) {
    CIS.METHODS <- c("flat", "raw")
    method <- pmatch(cis.method, CIS.METHODS)
    if (is.na(method)) {
      stop("invalid method")
    }
    if (method == -1) {
      stop("ambiguous method")
    }
    
    print("Working on cis data");
    ##
    ## Work on the cis chromosome 
    ##
    cis.data <- data[data[DATA.CHR] == chr,]
    cis.data.shuffle <- data[data[DATA.CHR] == chr,]
    cis.dom <- NULL
    use.raw <- TRUE
    if (CIS.METHODS[method] == "flat") {
      print("  This feature is not currently part of this implementation")
    }
    print("  Using the raw read counts")

    print(paste("Window size:", window.size))
              
    if (!is.null(mask.start) & !is.null(mask.end)) {
      ## we want "mask.start" to "round" to the
      ## nearest 5' fragment start.
      temp.start <- 0
      if (any(cis.data[DATA.POS] == mask.start)) {
        temp.start <- which(cis.data[DATA.POS] == mask.start)[[1]]
      } else {
        temp.start <- which(cis.data[DATA.POS] >= mask.start)[[1]]
        if (temp.start > 1) {
          temp.start <- temp.start - 1
        }
        new.mask.start <- cis.data[temp.start,DATA.POS]
        print("  mask was rounded to nearest 5' fragment start")
        print(paste("  ", mask.start, "is now", new.mask.start))
        mask.start <- new.mask.start
      }
      temp.end   <- max(which(cis.data[DATA.POS] <= mask.end))
      
      print("  Masking out...")
      if (invert.mask) {
        print(paste("... everything before (and including)", mask.start))
        print(paste("    and everything after (and including)", mask.end))
        print(paste("    on chromosome", chr))
        
        cis.data.shuffle <- cis.data[temp.start:temp.end,]
        cis.data <- cis.data[temp.start:temp.end,]
        
      } else {
        print(paste("  ...",
                    mask.start, "to",
                    mask.end, "on chromosome",
                    chr, "will be masked out"))
        
        cis.data.shuffle <- cis.data[c(1:temp.start, temp.end:nrow(cis.data)),]
        cis.data[temp.start:temp.end,TAB.READ.COUNT] <- 0
        
      }
    }
      
    cis.dom <- getDomains(cis.data, cis.data.shuffle, iterations=iterations,
                          window=window.size,
                          raw=use.raw,
                          fdr=fdr,
                          fdr.prob=fdr.prob,
                          fixed.cutoff=fixed.cutoff)
    ##if (print.results) {
    ##  print(cis.dom)
    ##}
  }

  if (cis.only) {
    if (!is.null(out.file)) {
      write.table(cis.dom, out.file, quote=F, col=F, row=F, sep="\t")
    }
    
    return(cis.dom)
  }
  
  if (!cis.only) {
    ## trans
    TRANS.METHODS <- c("flat", "raw")
    method <- pmatch(trans.method, TRANS.METHODS)
    if (is.na(method)) {
      stop("invalid method")
    }
    if (method == -1) {
      stop("ambiguous method")
    }

    print("Working on trans data")
    use.raw <- TRUE
    if (TRANS.METHODS[method] == "flat") {
      print("  This method is not part of this implementation")
    }
    print("  Using the raw read counts")
    
    trans.data <- data[data[DATA.CHR] != chr,]
    trans.dom <- NULL

    print(paste("Window size:", window.size))

    trans.dom <- getDomains(trans.data, trans.data, iterations=iterations,
                            window=window.size,
                            raw=use.raw,
                            fdr=fdr,
                            fdr.prob=fdr.prob,
                            fixed.cutoff=fixed.cutoff)
    ##if(print.results) {
    ##  print(trans.dom)
    ##}
      
  }

  if (trans.only) {
    if (!is.null(out.file)) {
      write.table(trans.dom, out.file, quote=F, col=F, row=F, sep="\t")
    }
    
    return(trans.dom)
  }
          
  ## package everthing up for printing to file.
  if (!is.null(out.file)) {
    if (!trans.only) {
      all.dom <- rbind(cis.dom, trans.dom)
      cis <- ifelse(all.dom[,1] == chr, 1, 0)
      output.df <- data.frame(all.dom, cis=cis)
      
      ## write the data to a file
      write.table(output.df, out.file, quote=F, col=F, row=F, sep="\t")
    }
  }
  
  return(list(cis = cis.dom, trans = trans.dom))
}

