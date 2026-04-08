 getPeakCRPF <- function(peakC.obj, fragmented.genome, re.length){
  #' PeakC does only include the position of the end of the fragment.
  #' Therefore, the r3c seq fragment table(fragmented genome) is used 
  #' to get corresponding start positions.
  #' Since peakC counts the re to the fragment end while r3c seq doesnt
  #' Cutter length has to be substracted from the postition of the peakC 
  #' fragment end
  #?? CAROLIN not all positions are found
  
  # New column containing the end positition with the cutter length sustracted 
  peakC.obj$dbR$frag_pos_minus_re <-peakC.obj$dbR$frag_pos - re.length
  
  # Extract positions of peaks 
  dbr.peaks <- peakC.obj$dbR[peakC.obj$dbR$frag_pos %in% peakC.obj$peak,]
  
  # Match peakC fragment ends with fragment ends from fragmented genome
  peaks.ends <- peakC.obj$peak[(peakC.obj$peak - re.length) %in% 
                                 end(fragmented.genome) ]-re.length
  
  #' Extract the full peakC fragment postions by matching peak ends 
  #' with the fragment ends of the fragmented genomes
  frag.genome.peaks <- fragmented.genome[end(fragmented.genome) %in% peaks.ends]
  
  
  # check, TODO for dev
  #dbr.peaks$frag_pos %>% length()
  #dbr.peaks$frag_pos %>% unique() %>% length()
  
  frag.genome.peaks.df <- as.data.frame(frag.genome.peaks) %>%
    mutate(tool = 'peakC') %>% 
    mutate(nReads = 0)  # TODO delete later
  
  # TODO unnötige schritte mit dbr peaks, da score nicht verwendet
  # ??Carolin brauchen wir peakC score
  # Name is not correct in meaning, but its necessary for merge
  frag.genome.peaks.df$frag_pos_minus_re <- frag.genome.peaks.df$end
  peakC.obj.full.df <- merge(dbr.peaks, frag.genome.peaks.df, by = 'frag_pos_minus_re') %>% 
    select(seqnames, start, end, nReads, tool)
  
  return(peakC.obj.full.df)
}

getUnion <- function(intervalls = list()){
  # Create GRanges Object
  # Apply schreibt nicht in global env
  intervalls.gr <- list()
  
  for ( i in seq_along(intervalls)){
    intervalls.gr[[i]] <- makeGRangesFromDataFrame(intervalls[[i]], keep.extra.columns = T)
  }
  # Give names to the list
  names(intervalls.gr) <- names(intervalls)
  
  # Union gr objects
  union.gr <- Reduce(GenomicRanges::union, intervalls.gr)
  
  # Look into data
  union.gr %>% head()
  
  # Plot to check if it worked out
  # Prepare data
  union.df <- union.gr %>% 
    as.data.frame() %>% 
    select(seqnames, start, end) %>%
    mutate(nReads = 0) %>% 
    mutate(tool = 'union')
  
  
  # Merge it together with intervalls
  union.merged <- rbind(Reduce(rbind, intervalls), union.df)
  
  # Bring union to bottom in plot by sorting levels
  union.merged$tool <- fct_relevel(union.merged$tool, "union", after = Inf)
  
  return(list(union.gr = union.gr,
              union.df = union.merged))
  # EOF
}

getUnion.grl <- function(intervalls = GRangesList()){
  # Get union of list elements
  union.gr <- Reduce(GenomicRanges::union, intervalls)
  union.gr$nReads <- 0
  union.gr$tool <- 'union'
  # Look into data
  union.gr %>% head()
  return(union.gr)
}



getIntersect <- function(intervalls = list()){
  # Find overlap
  # Create GRanges Object
  # Apply schreibt nicht in global env
  intervalls.gr <- list()
  for ( i in seq_along(intervalls)){
    intervalls.gr[[i]] <- makeGRangesFromDataFrame(intervalls[[i]], keep.extra.columns = T)
  }
  
  names(intervalls.gr) <- names(intervalls)
  
  # Find overlap
  intersect <- Reduce(GenomicRanges::intersect, intervalls.gr) %>%
    as.data.frame() %>%
    select(seqnames, start, end) %>%
    mutate(nReads = 0) %>%
    mutate(tool = 'intersect total')  
  
  if (dim(intersect)[1] == 0)  warning('No intersect found')  
  
  # Merge overlap with individual intervalls
  # overlap$tool <- factor(overlap$tool, levels = levels.ordered)
  intersect.merged <- rbind(Reduce(rbind, intervalls), intersect)
  
  # Order Factors overlap total has to be the last
  library(forcats)
  
  # TODO add this in plot stiff?? or delete
  intersect.merged$tool <- fct_relevel(intersect.merged$tool, "intersect total", after = Inf)
  
  levels(intersect.merged$tool) <- c(levels(intersect.merged$tool),  "intersect total")
  #levels.ordered <- fct_relevel(overlap.merged$tool %>% levels, "overlap total", after = Inf)
  
  
  return(union.df = intersect.merged)
}

getIntersect.grl <- function(intervalls = GRangesList(), name = NULL, out = FALSE){
  # Get intersect
  # intersect <- data.frame(seqnames = character(), start = integer(), 
  #                                    end = integer(),
  #                                    nReads = integer(), tool = factor())
  intersect <- Reduce(GenomicRanges::intersect, intervalls) #%>%

  
  if (length(intersect) == 0) {
    warning('No intersect found')  
    gr <-  makeGRangesFromDataFrame(data.frame(seqnames='chr1',start=0,end=0, nReads = -1, tool = ifelse(is.null(name), 'intersect', name)), keep.extra.columns = T)
    return(gr)
  } else if (length(intersect) == 1  && intersect$nReads == -1){
    intersect$tool <- ifelse(is.null(name), 'intersect', name)
    return(intersect)
    
   } else {
    intersect$nReads <- 0
    intersect$significance <- 1  # binary -> could be not binary
    intersect$tool <- ifelse(is.null(name), 'intersect', name)
    return(intersect)
  }
}


getCustomizedPairwiseIntersect <- function(intervalls = list(), intersect.of.interest = character()){
  intervalls.gr <- list()
  for ( i in seq_along(intervalls)){
    #assign(paste0(names(intervalls)[i], '.gr'), makeGRangesListFromDataFrame(intervalls[[i]]))
    intervalls.gr[[i]] <- makeGRangesFromDataFrame(intervalls[[i]], keep.extra.columns = T)
  }
  
  names(intervalls.gr) <- names(intervalls)
  
  
  # Pairwise overall of interest
  # Plot selected paiwise contrast
  i.o.i <-  strsplit(x = intersect.of.interest, split = '\\+')[[1]]
  i.o.i.gr <- GenomicRanges::intersect(intervalls.gr[[i.o.i[1]]], 
                                       intervalls.gr[[i.o.i[2]]])
  
  
  i.o.i.df <- i.o.i.gr %>%
    as.data.frame() %>% 
    select(seqnames, start, end) %>%
    mutate(tool = intersect.of.interest) %>%
    mutate(nReads = 0)
  
  i.o.i.merged <- rbind(Reduce(rbind, intervalls), i.o.i.df)
  
  # Prepare data for plotting
  # Transform tools to factor
  levels <- i.o.i.merged$tool %>% as.factor() %>% levels()
  
  # Order the tools so that intersects are the last ones
  levels.ordered <- levels[order(str_count(levels, '\\+'))]
  
  # Merge overlap with individual intervalls
  i.o.i.merged$tool <- factor(i.o.i.merged$tool, levels = levels.ordered)
  return(i.o.i.merged)
}
getMajorityVote <- function(intervalls = list(), greater.than = int()){
  # Collect intervalls gr in list
  intervalls.gr <- list()
  
  # Transform intervalls to gr
  for ( i in seq_along(intervalls)){
    intervalls.gr[[i]] <- makeGRangesFromDataFrame(intervalls[[i]], keep.extra.columns = T)
  }
  
  # Give names to the list
  names(intervalls.gr) <- names(intervalls)
  
  
  concat.gr <- Reduce(rbind, intervalls) %>% makeGRangesFromDataFrame(., keep.extra.columns = T)
  cov <- coverage(concat.gr)
  
  mv.df <- as(cov, "GRanges") %>% 
    as.data.frame() %>% 
    mutate(nReads = 0) %>% 
    mutate(tool = paste0('majority vote >= ', as.character(greater.than))) %>%
    filter(score >= greater.than) %>%
    dplyr::select(seqnames, start, end, nReads, tool)
  
  
   mv.merged <- rbind(Reduce(rbind, intervalls), mv.df) # was turned back earlier
  
  # Bring union to bottom in plot by sorting levels
  mv.merged$tool <- fct_relevel(mv.merged$tool, paste0('majority vote >= ', as.character(greater.than)), after = Inf)
  return(list(mv = as(cov, "GRanges"),
              mv.df = mv.df))
}

getMajorityVote.grl <- function(intervalls = list(), greater.than = int()){
  # Collect intervalls gr in list
  intervalls.gr <- intervalls@unlistData
  
  # Transform intervalls to gr
  # for ( i in seq_along(intervalls)){
  #   intervalls.gr[[i]] <- makeGRangesFromDataFrame(intervalls[[i]], keep.extra.columns = T)
  # }
  
  # Give names to the list
  
  cov <- coverage(intervalls.gr)
  
  # TODO change to GRanges
  mv <- as(cov, "GRanges") #%>% 
  mv <- mv[mv$score >= greater.than,]
  
  if (length(mv) > 0){
    mv$tool <- paste0('majority vote >= ', as.character(greater.than))
    #mv <- mv[, 'tool']
    mv %>% as.data.frame() %>%
    mutate(nReads = 0) %>%
    mutate(tool = paste0('majority vote >= ', as.character(greater.than))) %>%
  #  filter(score >= greater.than) %>%
    dplyr::select(seqnames, start, end, nReads, tool)
    
    mv$significance <- 1  # binary
    mv <- makeGRangesFromDataFrame(mv, keep.extra.columns = T)
    
  } else {
    mv <- makeGRangesFromDataFrame(
      data.frame(seqnames='chr1', start=0, end=0, nReads = -1, significance = 0, 
                                              tool = paste0('majority vote >= ',
                                                            as.character(
                                                              greater.than))), 
      keep.extra.columns = TRUE)
  }
  return(mv)
}









# getMajorityVote.grl <- function(intervalls = list(), greater.than = int()){
#   # Collect intervalls gr in list
#   intervalls.gr <- intervalls@unlistData
#   
#   # Transform intervalls to gr
#   # for ( i in seq_along(intervalls)){
#   #   intervalls.gr[[i]] <- makeGRangesFromDataFrame(intervalls[[i]], keep.extra.columns = T)
#   # }
#   
#   # Give names to the list
#   
#   
#   cov <- coverage(intervalls.gr)
#   
#   # TODO change to GRanges
#   mv <- as(cov, "GRanges") #%>% 
#   mv <- mv[mv$score >= greater.than,]
#   
#   if (length(mv) > 0){
#     mv$tool <- paste0('majority vote >= ', as.character(greater.than))
#     mv <- mv[, 'tool']
#     # as.data.frame() %>% 
#     # mutate(nReads = 0) %>% 
#     # mutate(tool = paste0('majority vote >= ', as.character(greater.than))) %>%
#     # filter(score >= greater.than) %>%
#     # select(seqnames, start, end, nReads, tool) %>%
#     # makeGRangesFromDataFrame()
#     
#     return(mv)
#   }
# }
# 
# 
# getMajorityVote.grl <- function(intervalls = list(), greater.than = int()){
#   # Collect intervalls gr in list
#   intervalls.gr <- intervalls@unlistData
#   
#   # Transform intervalls to gr
#   # for ( i in seq_along(intervalls)){
#   #   intervalls.gr[[i]] <- makeGRangesFromDataFrame(intervalls[[i]], keep.extra.columns = T)
#   # }
#   
#   # Give names to the list
#   
#   
#   cov <- coverage(intervalls.gr)
#   
#   # TODO change to GRanges
#   mv <- as(cov, "GRanges") #%>% 
#   mv <- mv[mv$score >= greater.than,]
#   
#   if (length(mv) > 0){
#     mv$tool <- paste0('majority vote >= ', as.character(greater.than))
#     mv <- mv[, 'tool']
#     # as.data.frame() %>% 
#     # mutate(nReads = 0) %>% 
#     # mutate(tool = paste0('majority vote >= ', as.character(greater.than))) %>%
#     # filter(score >= greater.than) %>%
#     # select(seqnames, start, end, nReads, tool) %>%
#     # makeGRangesFromDataFrame()
#     
#     return(mv)
#   }
# }


plotConsens <- function(consens.df = data.frame(),
                        read.counts.gr = GRanges(),
                        xlim = c(NULL, NULL),
                        ymax = NULL,
                        colors = NULL){
  
  if (length(read.counts.gr) == 0){
    consens.df <- consens.df %>% mutate(y.pos = -as.numeric(tool))
    
    p <- ggplot() +
      geom_rect(data = consens.df, 
                aes(xmin = start, xmax = end, ymin =  y.pos-1, 
                    ymax = y.pos, fill = tool)) +
      coord_cartesian(xlim = xlim) +
      labs(x = 'Position') 
    return(p)
    
  } else {
    # format read count table   
    read.counts.df <- read.counts.gr %>% 
      as.data.frame() %>%
      select(-width, -strand)
    # TODO code soft
    # TODO melt
    colnames(read.counts.df)[ colnames(read.counts.df) =='reads'] <- paste0('read_counts') 
    
    
    if (!is.numeric(ymax)){
      line.height <- read.counts.df[4] %>% max()
      consens.df <- consens.df %>% 
        mutate(y.pos = -as.numeric(tool)*line.height*0.1+line.height*0.1)
      
      if ('condition' %in% colnames(read.counts.df)){
        p <- ggplot() +
          geom_rect(data = consens.df, 
                    aes(xmin = start, xmax = end, ymin =  y.pos-(line.height*0.1), 
                        ymax = y.pos, fill = tool)) +
          geom_line(data = read.counts.df,
                    aes(x = start, y = read_counts, color = condition)) ####
        return(p)
      } else {
        p <- ggplot() +
          geom_rect(data = consens.df, 
                    aes(xmin = start, xmax = end, ymin =  y.pos-(line.height*0.1), 
                        ymax = y.pos, fill = tool)) +
          geom_line(data = read.counts.df,
                    aes(x = start, y = read_counts)) ####
      }
      return(p)
      
      
    } else {
      line.height <- ymax
      consens.df <- consens.df %>% 
        mutate(y.pos = -as.numeric(tool)*line.height*0.1+line.height*0.1)
      
      
      if ('condition' %in% colnames(read.counts.df)){
        p <- ggplot() +
          geom_rect(data = consens.df, 
                    aes(xmin = start, xmax = end, ymin =  y.pos-(line.height*0.1), 
                        ymax = y.pos, fill = tool)) +
          geom_line(data = read.counts.df,
                    aes(x = start, y = read_counts, color = condition)) +
          coord_cartesian(xlim = xlim, ylim = c(min(consens.df$y.pos)-(line.height*0.1), ymax)) 
      } else {
        p <- ggplot() +
          geom_rect(data = consens.df, 
                    aes(xmin = start, xmax = end, ymin =  y.pos-(line.height*0.1), 
                        ymax = y.pos, fill = tool)) +
          geom_line(data = read.counts.df,
                    aes(x = start, y = read_counts)) +
          coord_cartesian(xlim = xlim, ylim = c(min(consens.df$y.pos)-(line.height*0.1), ymax)) 
      }
      
      if (!is.null(colors)){
        p <- p + scale_fill_manual(values = colors)
        
      }
      
      return(p)
    }
  }
}



plotConsens.grl <- function(consens.df = GRangesList(),
                            read.counts.gr = GRanges(),
                            xlim = c(NULL, NULL),
                            ymax = NULL,
                            colors = 'standard', void = T){
  
  # Transfer GRlist into dataframe
  consens.df = GenomicRanges::as.data.frame(consens.df@unlistData, row.names = NULL) %>%
    dplyr::select(any_of(c('seqnames', 'start', 'end', 'nReads', 'tool', 'significance'))) %>%
    mutate(tool = as.factor(tool))
  
  consens.df$tool <- factor(consens.df$tool, levels = unique(droplevels(consens.df$tool)))
  consens.df <- consens.df %>% mutate(y.pos = -as.numeric(tool))
  consens.df$significance <- as.factor(consens.df$significance)
  # consens.df$significance <- as.factor(consens.df$significance)
  
  
  # If no read counts should be plotted (grl)
  if (length(read.counts.gr) == 0){
    y_labels <- unique(consens.df[, c('tool', 'y.pos')])
    y_label_map <- setNames(as.character(y_labels$tool), y_labels$y.pos)
    
   if ('significance' %in% colnames(consens.df)){
     p <- ggplot(data = consens.df)
     
     # Conditional layer for 'known' tool
     if ('known' %in% levels(consens.df$tool)) {
       known_layer <- geom_rect(
         data = consens.df[consens.df$tool == 'known', ],
         aes(xmin = start, xmax = end, ymin = max(consens.df$y.pos),
             ymax = min(consens.df$y.pos) - 1), 
         fill = 'grey70', 
         show.legend = FALSE
       )
       p <- p + known_layer
     }
     
     # Add the core geom_rect for all data
     p <- p + geom_rect(
       aes(xmin = start, xmax = end, ymin = y.pos - 1,
           ymax = y.pos, fill = tool, alpha = significance),
       show.legend = FALSE
     )
      p <- p + 
     
     
       scale_alpha_manual(values = c(`0` = 0, `1` = 0.33, `2` = 0.66, `3` = 1)) +
       #scale_alpha(range = c(0.25, 1)) +
       coord_cartesian(xlim = xlim) +
       labs(x = 'Position') +
       scale_y_continuous(breaks = unique(consens.df$y.pos - 0.5),
                          labels = function(x) y_label_map[as.character(x + 0.5)])

     #}
     
   } else {
     p <- ggplot() +
       geom_rect(data = consens.df,
                 aes(xmin = start, xmax = end, ymin =  y.pos-1,
                     ymax = y.pos, 
                     fill = tool), show.legend = F) +
       # Remove the geom_text function, as it is no longer needed
       coord_cartesian(xlim = xlim) +
       labs(x = 'Position', y = 'Tool') +
       scale_y_continuous(breaks = unique(consens.df$y.pos - 0.5),
                          labels = function(x) y_label_map[as.character(x + 0.5)])
   }
    

    
  } else {
    # format read count table   
    read.counts.df <- read.counts.gr
    colnames(read.counts.df)[ colnames(read.counts.df) =='reads'] <- paste0('read_counts')
    
    
    if (!is.numeric(ymax)){
      line.height <- read.counts.df$read_counts %>% max()
      consens.df <- consens.df %>%
        mutate(y.pos = -as.numeric(tool)*line.height*0.1+line.height*0.1)
      
      if ('condition' %in% colnames(read.counts.df)){
        p <- ggplot() +
          geom_rect(data = consens.df,
                    aes(xmin = start, xmax = end, ymin =  y.pos-(line.height*0.1),
                        ymax = y.pos, fill = tool), show.legend = F) +
          geom_line(data = read.counts.df,
                    aes(x = start, y = read_counts, color = condition)) +
        coord_cartesian(xlim = xlim)
        
      } else {
        p <- ggplot() +
          geom_rect(data = consens.df,
                    aes(xmin = start, xmax = end, ymin =  y.pos-(line.height*0.1),
                        ymax = y.pos, fill = tool)) +
          geom_line(data = read.counts.df,
                    aes(x = start, y = read_counts)) +
          coord_cartesian(xlim = xlim, ylim = c(min(consens.df$y.pos)-(line.height*0.1), ymax))
        
      }
      #return(p)
      
      
    } else {
      line.height <- ymax
      consens.df <- consens.df %>%
        mutate(y.pos = -as.numeric(tool)*line.height*0.1+line.height*0.1)
      
      
      if ('condition' %in% colnames(read.counts.df)){
        p <- ggplot() +
          geom_rect(data = consens.df,
                    aes(xmin = start, xmax = end, ymin =  y.pos-(line.height*0.1),
                        ymax = y.pos, fill = tool)) +
          geom_line(data = read.counts.df,
                    aes(x = start, y = read_counts, color = condition)) +
          coord_cartesian(xlim = xlim, ylim = c(min(consens.df$y.pos)-(line.height*0.1), ymax))
      } else {
        p <- ggplot() +
          geom_rect(data = consens.df,
                    aes(xmin = start, xmax = end, ymin =  y.pos-(line.height*0.1),
                        ymax = y.pos, fill = tool)) +
          geom_line(data = read.counts.df,
                    aes(x = start, y = read_counts)) +
          coord_cartesian(xlim = xlim)
      }
      
      
      #return(p)
    }
  }
  
  if (!is.null(colors) & colors != 'standard'){
    p <- p + scale_fill_manual(values = colors)
    
  } else if (colors == 'standard'){
    levels <- unique((p$plot_env$consens.df$tool))

    levels.vec <- as.character(levels)
    names(levels.vec) <- levels.vec
    
    for (i in names(levels.vec)){
      if (is.na(i)){
      } else if (str_detect(i, '^known$')){
        levels.vec[i] <- 'black'
      } else if (str_detect(i, '^(single|rep).foursig.*')){
        levels.vec[i] <- '#3288bd'
      } else if (str_detect(i, '^(single|rep).peakc.*')){
        levels.vec[i] <- '#9c0142'
      } else if (str_detect(i, '^(single|rep).r3c.*')){
        levels.vec[i] <- '#fdae61'
      } else if (str_detect(i, '^(single|rep).r4cker.*')){
        levels.vec[i] <- 'forestgreen'
      } else if (str_detect(i, '^(single|rep).r4cker.*')){
        levels.vec[i] <- 'violet'
        # levels.vec[i] <- 'x'
      # } else if (str_detect(i, '^(inter|union)\\.(foursig|peakc)\\+(foursig|peakc).*')){
      #   levels.vec[i] <- '#5e4fa2'      
      # } else if (str_detect(i, '^(inter|union)\\.(foursig|r3c)\\+(foursig|r3c).*')){
      #   levels.vec[i] <- '#66c2a5'   
      # } else if (str_detect(i, '^(inter|union)\\.(peakc|r3c)\\+(peakc|r3c).*')){
      #   levels.vec[i] <- '#f46d43'      
        
      } else {
        levels.vec[i] <- 'grey20'      
      }
    }    

    
    p <- p + scale_fill_manual(values = levels.vec) +
      labs(x = paste0('Position ', consens.df$seqnames[1]),
          # y = '',
           color = 'Tool') +
      #theme(axis.text.y = element_blank()) +
      scale_y_continuous(breaks = unique(consens.df$y.pos - 0.5),
                         labels = function(x) y_label_map[as.character(x + 0.5)]) +
      theme(panel.grid = element_blank(), 
            panel.background = element_blank(),
            panel.grid.minor.y = element_line(color = 'grey', linetype = 1))
    # col = ifelse(str_detect(pattern = 'foursig',string = levels), 'blue',
    #        ifelse(str_detect(pattern = 'r3c', string = levels), 'orange', 'grey'))
    # p + scale_fill_manual(values = ifelse(grepl("foursig", levels), "blue",
    #          ifelse(grepl("peakc", levels), "red",
    #                 ifelse(grepl("r3c", levels), "yellow",
    #                      'grey'))))
    #
    # p
    
  }
  
  if (void){
    return(p)
  } else {
  return(list(p, consens.df))
  }
}




# plotConsens.grl <- function(consens.df = GRangesList(),
#                         read.counts.gr = GRanges(),
#                         xlim = c(NULL, NULL),
#                         ymax = NULL,
#                         colors = NULL){
#   
#   # If no read counts should be plotted (grl)
#   if (length(read.counts.gr) == 0){
#     
#     # Transfer GRlist into dataframe
#     consens.df = GenomicRanges::as.data.frame(consens.df@unlistData, row.names = NULL) %>%
#       select(seqnames, start, end, nReads, tool) %>%
#       mutate(tool = as.factor(tool))
#     consens.df <- consens.df %>% mutate(y.pos = -as.numeric(tool))
#     
#     p <- ggplot() +
#       geom_rect(data = consens.df, 
#                 aes(xmin = start, xmax = end, ymin =  y.pos-1, 
#                     ymax = y.pos, fill = tool)) +
#       coord_cartesian(xlim = xlim) +
#       labs(x = 'Position') 
#     return(p)
#     
#     
#   } else {
#     # format read count table   
#     read.counts.df <- read.counts.gr %>% 
#       as.data.frame() %>%
#       select(-width, -strand)
#     # TODO code soft
#     # TODO melt
#     colnames(read.counts.df)[ colnames(read.counts.df) =='reads'] <- paste0('read_counts') 
#     
#     
#     if (!is.numeric(ymax)){
#       line.height <- read.counts.df[4] %>% max()
#       consens.df <- consens.df %>% 
#         mutate(y.pos = -as.numeric(tool)*line.height*0.1+line.height*0.1)
#       
#       if ('condition' %in% colnames(read.counts.df)){
#         p <- ggplot() +
#           geom_rect(data = consens.df, 
#                     aes(xmin = start, xmax = end, ymin =  y.pos-(line.height*0.1), 
#                         ymax = y.pos, fill = tool)) +
#           geom_line(data = read.counts.df,
#                     aes(x = start, y = read_counts, color = condition)) ####
#         return(p)
#       } else {
#         p <- ggplot() +
#           geom_rect(data = consens.df, 
#                     aes(xmin = start, xmax = end, ymin =  y.pos-(line.height*0.1), 
#                         ymax = y.pos, fill = tool)) +
#           geom_line(data = read.counts.df,
#                     aes(x = start, y = read_counts)) ####
#       }
#       return(p)
#       
#       
#     } else {
#       line.height <- ymax
#       consens.df <- consens.df %>% 
#         mutate(y.pos = -as.numeric(tool)*line.height*0.1+line.height*0.1)
#       
#       
#       if ('condition' %in% colnames(read.counts.df)){
#         p <- ggplot() +
#           geom_rect(data = consens.df, 
#                     aes(xmin = start, xmax = end, ymin =  y.pos-(line.height*0.1), 
#                         ymax = y.pos, fill = tool)) +
#           geom_line(data = read.counts.df,
#                     aes(x = start, y = read_counts, color = condition)) +
#           coord_cartesian(xlim = xlim, ylim = c(min(consens.df$y.pos)-(line.height*0.1), ymax)) 
#       } else {
#         p <- ggplot() +
#           geom_rect(data = consens.df, 
#                     aes(xmin = start, xmax = end, ymin =  y.pos-(line.height*0.1), 
#                         ymax = y.pos, fill = tool)) +
#           geom_line(data = read.counts.df,
#                     aes(x = start, y = read_counts)) +
#           coord_cartesian(xlim = xlim, ylim = c(min(consens.df$y.pos)-(line.height*0.1), ymax)) 
#       }
#       
#       if (!is.null(colors)){
#         p <- p + scale_fill_manual(values = colors)
#         
#       }
#       
#       return(p)
#     }
#   }
# }
# 

getPrecisionRecallF1.gr.frag <- function(gold.standard, experiment, frag ){
  ov <- findOverlaps(frag, gold.standard)
  
  # Get the indices of the fragments hit by at least one peak
  fragment_hits <- subjectHits(ov)
  gold.standard.mapped <- gold.standard[fragment_hits,]
  
  true.pos <- sum(gold.standard.mapped %over%  experiment)
  precision <- true.pos / length(experiment)
  recall <- true.pos / length(gold.standard.mapped)
  f1 <- 2*((precision*recall)/(precision+recall))
  
  # Store results in vector
  precision.recall <- data.frame(precision = precision,
                                 recall = recall, 
                                 f1 = f1)
  return(precision.recall)
}
getPrecisionRecallF1Mcc.gr.frag <- function(gold.standard, experiment, frag) {
  ov <- findOverlaps(frag, gold.standard)
  
  # Get the indices of the fragments hit by at least one peak
  fragment_hits <- subjectHits(ov)
  gold.standard.mapped <- gold.standard[fragment_hits,]
  
  
  # True Positives
  true.pos <- as.numeric(sum(gold.standard.mapped %over% experiment))
  
  # False Positives
  false.pos <- as.numeric(length(experiment) - true.pos)
  
  # False Negatives
  false.neg <- as.numeric(sum(!gold.standard.mapped %over% experiment))
  
  # Precision, Recall, and F1 Score
  precision <- true.pos / length(experiment)
  recall <- true.pos / length(gold.standard.mapped)
  f1 <- 2 * ((precision * recall) / (precision + recall))
  
  # Matthews Correlation Coefficient
  # We assume all items in `frag` are relevant and mapped to evaluate their 'negative' correlations
  total_elements <- length(frag) # Total number of elements - TODO nearbait!!!!!!
  true.neg <- total_elements - (true.pos + false.pos + false.neg)
  
  numerator <- (true.pos * true.neg) - (false.pos * false.neg)
  denominator <- sqrt((true.pos + false.pos) * (true.pos + false.neg) * (true.neg + false.pos) * (true.neg + false.neg))
  
  if (denominator == 0) {
    mcc <- NA
  } else {
    mcc <- numerator / denominator
  }
  
  # Store results in data frame
  precision.recall.mcc <- data.frame(
    precision = precision,
    recall = recall,
    f1 = f1,
    mcc = mcc
  )
  return(precision.recall.mcc)
}

getPrecisionRecallF1MCC_AUCPR.gr.frag <- function(gold.standard, experiment, frag, scale_rate = FALSE, auc_cutoff = 0.5){
  if (length(experiment) == 1 && start(experiment)[1] == 0 && end(experiment)[1] == 0){
    message('Tool was not calling anything before. Therefore, dummy was written. Returning all values 0.')
    precision.recall.mcc.aucpr <- data.frame(
      precision = 0,
      recall = 0,
      f1 = 0,
      mcc = 0,
      auc_pr = 0
    )
  } else {
    # Sort out fragments outside of frag area 
    # TODO map this in this step?
    #experiment <- experiment[experiment %in% frag,]
    
    
    
  library(PRROC)
    # # Map experiemnt
    ov <- findOverlaps(frag, experiment)
    
    fragment_hits <- subjectHits(ov)
    
    experiment.mapped <- frag
    experiment.mapped$rate_total <- 0
    
    
    if (scale_rate != T){
      mcols(experiment.mapped[queryHits(ov),])$rate_total <- experiment[subjectHits(ov),]$rate_total
    } else {
      mcols(experiment.mapped[queryHits(ov),])$rate_total <- (experiment[subjectHits(ov),]$significance)#/3  # warum /3 -> passt dnn ja nicht mit den Schranken später
    }
    
  
    
  ov <- findOverlaps(frag, gold.standard)
  
  # Get the indices of the fragments hit by at least one peak
  fragment_hits <- unique(queryHits(ov))
  gold.standard.mapped <- frag[fragment_hits,]
  
  
  # True Positives
  true.pos <- as.numeric(sum(gold.standard.mapped %over% experiment.mapped[experiment.mapped$rate_total > auc_cutoff]))
  
  # False Positives
  false.pos <- as.numeric(length(experiment.mapped[experiment.mapped$rate_total > auc_cutoff] ) - true.pos)
  
  # False Negatives
  false.neg <- as.numeric(sum(!gold.standard.mapped %over% experiment.mapped[experiment.mapped$rate_total > auc_cutoff] ))
  
  # Precision, Recall, and F1 Score
  precision <- true.pos / length(experiment.mapped[experiment.mapped$rate_total > auc_cutoff] )
  recall <- true.pos / length(gold.standard.mapped)
  f1 <- 2 * ((precision * recall) / (precision + recall))
  
  # Matthews Correlation Coefficient
  total_relevant_regions <- length(unique(c(queryHits(ov), subjectHits(findOverlaps(frag, experiment.mapped)))))
  true.neg <- max(0, total_relevant_regions - (true.pos + false.pos + false.neg))
  
  numerator <- (true.pos * true.neg) - (false.pos * false.neg)
  denominator <- sqrt(max(1e-8, (true.pos + false.pos) * (true.pos + false.neg) * (true.neg + false.pos) * (true.neg + false.neg)))
  
  
  
  if (denominator == 0) {
    mcc <- NA
  } else {
    mcc <- numerator / denominator
  }
  
  
  # AUC
  # Map vectors for AUC
  # Get the indices of the fragments hit by at least one peak
  # Map gold.standard to fragments
  ov <- findOverlaps(frag, gold.standard)
  
  fragment_hits <- subjectHits(ov)
  
  gold.standard.mapped.frag <- frag
  gold.standard.mapped.frag$rate_total <- 0
  
  # set gt to 1
  mcols(gold.standard.mapped.frag[queryHits(ov),])$rate_total <- 1
  
  rm(ov, fragment_hits)
  
  
  
  summary(experiment.mapped$rate_total)
  summary(gold.standard.mapped.frag$rate_total)
  
  pos.scores <- experiment.mapped$rate_total[gold.standard.mapped.frag$rate_total == 1]  
  neg.scores <- experiment.mapped$rate_total[gold.standard.mapped.frag$rate_total == 0]  

  
  
  # TODO if this is not working set AUC PR to 0
  # prroc_curve <-pr.curve(scores.class0 = experiment$rate_total, weights.class0 = gold.standard.mapped.frag$rate_total, curve = TRUE)
  #try(auc_pr <- prroc_curve$auc.integral)
  result <- tryCatch({
    prroc_curve <-pr.curve(scores.class0 = pos.scores, scores.class1 = neg.scores, curve = TRUE)
    auc_pr <- prroc_curve$auc.integral
    plot(prroc_curve)
    list(success = TRUE, auc = auc_pr)
  }, error = function(e) {
    message("Ein Fehler ist aufgetreten: ", e$message)
    list(success = FALSE, auc = NA)
  })
  
  # Überprüfe das Ergebnis
  if (result$success) {
    print(paste("AUC Integral: ", result$auc))
  } else {
    print("AUC Integral konnte nicht berechnet werden.")
    auc_pr <- NA
  }
  
  
# Print the AUC-PR result
  #try(cat("AUC-PR:", auc_pr, "\n"))
  
  ###
  
 
  try(precision.recall.mcc.aucpr <- data.frame(
    precision = precision,
    recall = recall,
    f1 = f1,
    mcc = mcc,
    auc_pr = auc_pr
  ))
  
  }
  return(precision.recall.mcc.aucpr)
  
  #
  # ov <- findOverlaps(frag, gold.standard)
  # fragment_hits <- queryHits(ov)
  # 
  # gold.standard.mapped <- frag
  # gold.standard.mapped$significance <- 0
  # mcols(gold.standard.mapped[fragment_hits,])$significance <- 1
  # 
  # rm(ov)
  # 
  # # Map experiments to fragments
  # ov <- findOverlaps(frag, experiment)
  # fragment_hits <- queryHits(ov)
  # 
  # experiment.mapped <- frag
  # # Retrieve actual prediction scores from `experiment`
  # experiment.mapped$significance <- 0
  # mcols(experiment.mapped[fragment_hits,])$significance <- experiment[subjectHits(ov),]$significance
  # predicted_binary <- ifelse(experiment.mapped$significance >= 1, 1, 0)
  # 
  # # Create vectors
  # ground.truth <- gold.standard.mapped$significance  # True binary labels
  # #experiment.scores <- experiment.mapped$significance  # Predicted scores
  # 
  # # Validate dimension and score range consistency
  # 
  # # Generate the precision-recall curve with PRROC
  # pr_curve <- pr.curve(
  #   scores.class1 = experiment.scores[ground.truth == 1],  # Scores for actual positives
  #   scores.class0 = experiment.scores[ground.truth == 0],  # Scores for actual negatives
  #   curve = TRUE
  # )
  # 
  # # Extract AUCPR value
  # auc_pr <- pr_curve$auc.integral
  # print(paste("The AUCPR is:", auc_pr))
  # 
  # # Plot the Precision-Recall curve
  # plot(pr_curve, main = "Precision-Recall Curve with AUCPR")
  # 
  
}


# Dev --------------------------------------------------------------------------
# # Create individual GRanges objects
# # dev ####
# gr1 <- GRanges(seqnames = "chr1",
#                ranges = IRanges(start = seq(0, 30, 5), end = (seq(0, 30, 5)+1)),
#                significance = c(1, 2, 3, 0, 0, 0, 0), tool = 'a')
# 
# gr2 <- GRanges(seqnames = "chr1",
#                ranges = IRanges(start = seq(0, 30, 5), end = (seq(0, 30, 5)+1)),
#                significance = c(0, 2, 3, 1, 0, 0, 0), tool = 'b')
# 
# gr3 <- GRanges(seqnames = "chr1",
#                ranges = IRanges(start = seq(0, 30, 5), end = (seq(0, 30, 5)+1)),
#                significance = c(0, 0, 0, 0, 1, 1, 1), tool = 'c')
# 
# gr4 <- GRanges(seqnames = "chr1",
#                ranges = IRanges(start = seq(0, 30, 5), end = (seq(0, 30, 5)+1)),
#                significance = c(0, 2, 3, 0, 0, 2, 1), tool = 'd')
# 
# frag <- GRanges(seqnames = "chr1",
#                 ranges = IRanges(start = seq(0, 30, 5),
#                                  end = (seq(0, 30, 5)+1)))
# 
# 
# 
# # Combine them into a GRangesList
# intervalls <- GRangesList(a = gr1, b = gr2, c = gr3, d = gr4)
# intervalls
# plotConsens.grl(intervalls)
# 
# # Creates weights object
# weights <- data.frame(weights = c(1,2,3,4))
# rownames(weights) <- letters[1:4]
# 
# 
# 
# # Add weights and rating to dataframe
# # (weights muss nicht sein, ohne effizienter, aber erstmal zur kontrolle mitnehmen)
# for (n in rownames(weights)){
#   intervalls[[n]]$weight <- weights[n, ]
#   intervalls[[n]]$rating <- intervalls[[n]]$significance*intervalls[[n]]$weight
# }
# 
# intervalls
# 
# 
# 
# # Calculate together
# for (inter in seq_along(intervalls)){
#   print(intervalls[[inter]])
#   mcols(frag)
#   mcols(frag)[[paste0('rate_', names(intervalls[inter]))]] <- intervalls[[inter]]$rating
#   
#   }
# rate.sum <- mcols(frag) %>% as.data.frame() %>% apply(., 1, sum)
# 
# 
# 
# maj <- frag[rate.sum > 10,]
# maj$significance <- 3
# maj$tool <- 'maj'
# intervalls[['maj']] <- maj
# 
# plotConsens.grl(intervalls)
# End dev ####



getMajorityVote.grl.sig <- function(intervalls = GRangesList(),
                                   weights = character(),
                                   greater.than = 10,
                                   frag = GRanges(),
                                   out = FALSE,
                                   mode = 'base'){
  
  # Extract intervalls included in weights
  intervalls <- intervalls[names(intervalls) %in%
                             paste0(rownames(weights), '.condition')]  # TODO not only condition
  # TODO abfangen, wenn ein nicht im anderen 
  
  # Add weights and rating to intervalls
  for (i in 1: nrow(weights)){
    print(paste0('weight: ', rownames(weights)[i], ' weights'))
    gr <- intervalls[startsWith(prefix =  paste0(rownames(weights)[i], '.condition'), x = names(intervalls))][[1]]
    mcols(intervalls[startsWith(prefix =  paste0(rownames(weights)[i], '.condition'), x = names(intervalls))][[1]])$weight <- weights[i,]  # TODO vlt. besser über name
    
    print(weights[i,] )
    #print(gr$significance[1])
    # Calculate rating
    mcols(intervalls[startsWith(prefix =  paste0(rownames(weights)[i], '.condition'), x = names(intervalls))][[1]])$rating <- gr$significance*weights[i,] 
  }
  
  # Map on fragments
  for (i in rownames(weights)){
    tmp <- paste0(i, '.condition')
    
    # Add significance information to virtual fragment library (vfl)
    # Map fragments onto vfl
    ov <- findOverlaps(frag, intervalls[[tmp]])
    
    # Get the indices of the fragments hit by at least one peak
    fragment_hits <- subjectHits(ov)
    
    # Extract the fragments hit
    mcols(frag)$significance <- 0
    mcols(frag)$rating <- 0
    
    frag[queryHits(ov)]$significance <-  mcols(intervalls[[tmp]])$significance[subjectHits(ov)]
    frag[queryHits(ov)]$rating <-  mcols(intervalls[[tmp]])$rating[subjectHits(ov)]
    # neuie mcol definieren
    
    frag$tool <- intervalls[[tmp]]$tool[1]
    intervalls[[tmp]] <- frag
    print(paste0('Map on frags done for ', tmp, '.'))
  }
  
  # warum zweites frag? -> damit ich da nicht zu viele metadaten reinhaue, aber warum
  
  # seit sig to 0 since only rating is used afterwards
  frag$significance <- 0
  # Calculate together
  for (inter in seq_along(intervalls)){
    
    #mcols(frag)[[paste0('rate_', names(intervalls[inter]))]] <- 'x'
    mcols(frag)[[paste0('rate_', names(intervalls[inter]))]] <- intervalls[[inter]]$rating
    
    # mcols(frag)[[paste0('rate_', names(intervalls[inter]))]] <- intervalls[[inter]]$rating
  }
  rate.sum <- mcols(frag) %>% 
    as.data.frame() %>%
    dplyr::select(starts_with('rate')) %>%
    apply(1, function(x) sum(unlist(x)))

  if (out == FALSE){
    maj <- frag[rate.sum > greater.than,]
    mcols(maj)$rate_total <- rate.sum[rate.sum > greater.than]
  } else {
    maj <- frag
    mcols(maj)$rate_total <- rate.sum
  }
  
  if (length(maj) > 0 && sum(maj$rate_total > greater.than) > 0){
    mcols(maj[maj$rate_total > greater.than,])$significance <- 3 # output in reality binary
  } else if (mode == 'AUC'){
    maj <- maj
  } 
  
  else {
    maj <- dummy.gr
  }
  # TODO
  maj$tool <- 'maj'
  # intervalls[['maj']] <- maj
  return(maj)
}

# Function to read a FASTQ file and determine the read length
# Function to read a FASTQ file and determine the read length
determine_read_length_from_fastq <- function(fastq_file) {
  # Open the FASTQ file for reading
  con <- file(fastq_file, "r")
  # Read the first sequence line  
  read_length <- NA
  while(TRUE) {
    # Read and discard the header line
    header <- readLines(con, n = 1)
    if(length(header) == 0) break # End of file
    
    # Read the sequence line
    sequence <- readLines(con, n = 1)
    # Determine length from the first sequence line
    read_length <- nchar(sequence)
    break # We only need the first sequence
    
    # Read and discard the plus line
    plus_line <- readLines(con, n = 1)
    # Read and discard the quality score line
    quality_score <- readLines(con, n = 1)
  }
  
  # Close the connection
  close(con)
  
  return(read_length)
}



