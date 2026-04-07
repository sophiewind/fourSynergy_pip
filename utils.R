readBedGraph <- function(ia) {
  bgs <- GRangesList()
  if (!is.null(ia@tracks)) {
    for (i in ia@metadata$conditionRep) {
      cond <- read.delim(
        paste0(
          ia@tracks, ia@metadata$condition, "_", i,
          "_sorted.bedGraph"
        ),
        header = FALSE
      ) %>%
        `colnames<-`(c("seqnames", "start", "end", "reads"))
      cond <- cond[startsWith(cond$seqnames, "chr"), ] %>%
        makeGRangesFromDataFrame(keep.extra.columns = TRUE)
      bgs[[paste0(ia@metadata$condition, "_", i)]] <- cond
    }

    for (i in ia@metadata$controlRep) {
      ctrl <- read.delim(
        paste0(
          ia@tracks, ia@metadata$control, "_", i,
          "_sorted.bedGraph"
        ),
        header = FALSE
      ) %>%
        `colnames<-`(c("seqnames", "start", "end", "reads"))
      ctrl <- ctrl[startsWith(ctrl$seqnames, "chr"), ] %>%
        makeGRangesFromDataFrame(keep.extra.columns = TRUE)
      bgs[[paste0(ia@metadata$control, "_", i)]] <- ctrl
    }
  }
  return(bgs)
}

trackPlots <- function(sia, tool = 'base', highlight_regions = NULL){
  bgs <- readBedGraph(sia)
  tmp <- sia@vfl

  # Collect reads
  collect <- list()
  for (i in seq(1, length(bgs))){
    ov <- findOverlaps(tmp, bgs[[i]])
    tmp$reads <- 0
    tmp[queryHits(ov)]$reads <- bgs[[i]][subjectHits(ov)]$reads
    collect[[names(bgs)[i]]] <- tmp$reads
  }
  coll.df <- as.data.frame(collect)

  # Average reads
  cond.reads <- coll.df %>%
    dplyr::select(., matches(sia@metadata$condition)) %>%
    mutate(mean_cond = rowMeans(.)) %>%
    dplyr::select(mean_cond)
  ctrl.reads <- coll.df %>%
    dplyr::select(., matches(sia@metadata$control)) %>%
    mutate(mean_ctrl = rowMeans(.)) %>%
    dplyr::select(mean_ctrl)

  tmp$cond_reads <- cond.reads
  tmp$ctrl_reads <- ctrl.reads


  if (!is.null(highlight_regions)){
      if (str_detect(highlight_regions,
                     "chr[1-9XYM]+\\:\\d+\\-\\d+(,\\s*chr[1-9XYM]+:\\d+-\\d+)*")){
          reg <- highlight_regions %>%
              stringr::str_split_1(',') %>%
              trimws() %>%
              as.data.frame() %>%
              separate('.', into = c("seqnames", "start", "end")) %>%
              makeGRangesListFromDataFrame()
      } else if (endsWith(highlight_regions, '.bed')){
          b.f <- read.delim(highlight_regions, header = FALSE)
          if (b.f[1,2] == "start"){
              b.f <- read.delim(highlight_regions, header = TRUE)
          }
          b.f <- b.f [,seq(1, 3)]
          colnames(b.f) <- c('seqnames', 'start', 'end')
          reg <- b.f %>%
              makeGRangesListFromDataFrame()
      } else {
          stop("The regions have to be either provided as string with following ",
               "pattern: 'chr3:33000000-34000000' or as comma separated string ",
               "(chr3:33000000-34000000, chr3:35000000-35500000) or as valid ",
               ".bed file.")
      }
  }

  tmp %>% head()

  if (tool == 'base'){
      # Import peaks
      for (p in c('rep.peakc_21.condition', 'rep.foursig_1.condition',
                  'rep.r3c_2000.condition', 'rep.r4cker_nearbait.condition',
                  'rep.peakc_21.control', 'rep.foursig_1.control',
                  'rep.r3c_2000.control', 'rep.r4cker_nearbait.control')){
          if (endsWith(p, 'condition')){
              ia <- sia@expInteractions[[p]]
          } else {
              ia <- sia@ctrlInteractions[[p]]
          }

          ov <- findOverlaps(tmp, ia)
          mcols(tmp)[[p]] <- 0
          mcols(tmp[queryHits(ov)])[[p]]<- ia[subjectHits(ov)]$significance
      }
      tmp.df <- as.data.frame(tmp, keep.extra.columns = TRUE)


      p.cond <-
          tmp.df %>%  dplyr::select(start, end, mean_cond, rep.r4cker_nearbait.condition,
                                    rep.peakc_21.condition, rep.foursig_1.condition,
                                    rep.r3c_2000.condition) %>%
          pivot_longer(!c('mean_cond', 'start', 'end'), names_to = "peak") %>%
          mutate(value = as.factor(value)) %>%
          ggplot()
      if (!is.null(highlight_regions)){
          p.cond <- p.cond + geom_rect(data = as.data.frame(reg),
                                       mapping = aes(xmin =  start, ymin =  0, xmax = end, ymax = Inf),
                                       alpha = 0.5)
      }
      p.cond <- p.cond + geom_segment(aes(x = start, xend = end, y = 0, yend = mean_cond), color = 'grey40') +
          ylim(c(0, 3000)) +
          geom_segment(aes(x = start, xend = end, y = 0, yend = mean_cond, alpha = value,
                           color = peak)) +
          scale_alpha_manual(values = c(0,1,1,1)) +
          # scale_color_manual(values = c('peak_r3c' = '#fdae61', 'peak_r4' = 'forestgreen',
          #                               'peak_peakc' = '#9c0142', 'peak_fs' = '#3288bd')) +
          facet_grid(.~peak) +
          scale_color_manual(values = c('rep.r3c_2000.condition' = '#fdae61', 'rep.r4cker_nearbait.condition' = 'forestgreen',
                                        'rep.peakc_21.condition' = '#9c0142', 'rep.foursig_1.condition' = '#3288bd')) +
          guides(alpha = "none", color = 'none') +
          labs(y = 'Counts', x = 'Position') +
          theme(panel.background = element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank()) +
          geom_vline(aes(xintercept = start(sia@vp)), linetype = 'dashed')


      p.ctrl <- tmp.df %>%  dplyr::select(start, end, mean_ctrl, rep.r4cker_nearbait.control,
                                          rep.peakc_21.control, rep.foursig_1.control,
                                          rep.r3c_2000.control) %>%
          pivot_longer(!c('mean_ctrl', 'start', 'end'), names_to = "peak") %>%
          mutate(value = as.factor(value)) %>%
          ggplot()
      if (!is.null(highlight_regions)){
          p.ctrl <- p.ctrl + geom_rect(data = as.data.frame(reg),
                                       mapping =
                                           aes(xmin =  start, ymin =  0, xmax = end, ymax = Inf),
                                       alpha = 0.5)
      }
      p.ctrl <- p.ctrl+ geom_segment(aes(x = start, xend = end, y = 0, yend = mean_ctrl), color = 'grey40') +
          ylim(c(0, 3000)) +
          geom_segment(aes(x = start, xend = end, y = 0, yend = mean_ctrl, alpha = value,
                           color = peak)) +
          scale_alpha_manual(values = c(0,1,1,1)) +
          # scale_color_manual(values = c('peak_r3c' = '#fdae61', 'peak_r4' = 'forestgreen',
          #                               'peak_peakc' = '#9c0142', 'peak_fs' = '#3288bd')) +
          facet_grid(.~peak) +
          scale_color_manual(values = c('rep.r3c_2000.control' = '#fdae61', 'rep.r4cker_nearbait.control' = 'forestgreen',
                                        'rep.peakc_21.control' = '#9c0142', 'rep.foursig_1.control' = '#3288bd')) +
          guides(alpha = "none", color = 'none') +
          labs(y = 'Counts', x = 'Position') +
          theme(panel.background = element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank()) +
          geom_vline(aes(xintercept = start(sia@vp)), linetype = 'dashed')
  } else {
      tmp$ens_cond <- sia@expConsensus$significance
      tmp$ens_ctrl <- sia@ctrlConsensus$significance
      tmp.df <- as.data.frame(tmp, keep.extra.columns = TRUE)

      p.cond <-
          tmp.df %>%  dplyr::select(start, end, mean_cond, ens_cond) %>%
          pivot_longer(!c('mean_cond', 'start', 'end'), names_to = "peak") %>%
          mutate(value = as.factor(value)) %>%
          ggplot()
      if (!is.null(highlight_regions)){
          p.cond <- p.cond + geom_rect(data = as.data.frame(reg),
                                       mapping =
                                           aes(xmin =  start, ymin =  0, xmax = end, ymax = Inf),
                                       alpha = 0.5)
      }
      p.cond <- p.cond + geom_segment(aes(x = start, xend = end, y = 0, yend = mean_cond), color = 'grey40') +
          ylim(c(0, 3000)) +
          geom_segment(aes(x = start, xend = end, y = 0, yend = mean_cond, alpha = value,
                           color = peak)) +
          scale_alpha_manual(values = c(0,1,1,1)) +
          # scale_color_manual(values = c('peak_r3c' = '#fdae61', 'peak_r4' = 'forestgreen',
          #                               'peak_peakc' = '#9c0142', 'peak_fs' = '#3288bd')) +
          facet_grid(.~peak) +
          scale_color_manual(values = c('ens_cond' = 'purple')) +
          guides(alpha = "none", color = 'none') +
          labs(y = 'Counts?', x = 'Position') +
          theme(panel.background = element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank()) +
          geom_vline(aes(xintercept = start(sia@vp)), linetype = 'dashed')


      p.ctrl <- tmp.df %>%  dplyr::select(start, end, mean_ctrl, ens_ctrl) %>%
          pivot_longer(!c('mean_ctrl', 'start', 'end'), names_to = "peak") %>%
          mutate(value = as.factor(value)) %>%
          ggplot()
      if (!is.null(highlight_regions)){
          p.ctrl <- p.ctrl + geom_rect(data = as.data.frame(reg),
                                       mapping =
                                           aes(xmin =  start, ymin =  0, xmax = end, ymax = Inf),
                                       alpha = 0.5)
      }
      p.ctrl <- p.ctrl+ geom_segment(aes(x = start, xend = end, y = 0, yend = mean_ctrl), color = 'grey40') +
          ylim(c(0, 3000)) +
          geom_segment(aes(x = start, xend = end, y = 0, yend = mean_ctrl, alpha = value,
                           color = peak)) +
          scale_alpha_manual(values = c(0,1,1,1)) +
          # scale_color_manual(values = c('peak_r3c' = '#fdae61', 'peak_r4' = 'forestgreen',
          #                               'peak_peakc' = '#9c0142', 'peak_fs' = '#3288bd')) +
          facet_grid(.~peak) +
          scale_color_manual(values = c('ens_ctrl' = 'purple')) +
          guides(alpha = "none", color = 'none') +
          labs(y = 'Counts?', x = 'Position') +
          theme(panel.background = element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank()) +
          geom_vline(aes(xintercept = start(sia@vp)), linetype = 'dashed')
  }
  return(list(p.cond, p.ctrl))
}

get_gene_selection <- function(id, sia) {
    TxDb <- switch(sia@metadata$organism,
                   "mm9" = TxDb.Mmusculus.UCSC.mm9.knownGene,
                   "mm10" = TxDb.Mmusculus.UCSC.mm10.knownGene,
                   "hg19" = TxDb.Hsapiens.UCSC.hg19.knownGene,
                   "hg38" = TxDb.Hsapiens.UCSC.hg38.knownGene)

    # Provide gene selection
    genes <- genes(TxDb)
    area <- GRanges(seqnames = seqnames(sia@vp)[1],
                    IRanges(start(sia@vfl)[1], end(sia@vfl)[length(sia@vfl)]))
    overlapping_genes <- findOverlaps(genes, area)
    genes.area <- genes[queryHits(overlapping_genes)]
    if (startsWith(sia@metadata$organism, 'mm')){
        org <- org.Mm.eg.db
    } else {
        org <- org.Hs.eg.db
    }
    trans <- clusterProfiler::bitr(geneID = genes.area$gene_id,
    fromType = "ENTREZID", "SYMBOL", org)
  renderUI({selectizeInput(paste0("genes_select_", id), "Genes:",
                           choices = trans$SYMBOL,
                           multiple = TRUE)
  })
}


#' Title
#'
#' @param object
#' @param value
#'
#' @returns
#' @export
#'
#' @examples
setDifferential <- function(object, value) {
  if (inherits(value, "DESeqResults")) {
    object@differential <- value
  } else {
    stop("Value muss ein DESeqResults-Objekt sein")
  }
  return(object)
}

#' Title
#'
#' @param object
#' @param value
#'
#' @returns
#' @export
#'
#' @examples
setDds <- function(object, value) {
    if (inherits(value, "DESeqDataSet")) {
        object@dds <- value
    } else {
        stop("Value must be DESeqDataSet object.")
    }
    return(object)
}

#info_file <- "/host/Datasets/Geeven_sox/info.yaml"
check_info <- function(info_file) {
    if (is.null(info_file) || !(grepl("\\.yaml$", info_file, ignore.case = TRUE))) {
        stop("Info file is missing or not of type .yaml")
    }
    pattern <- "^[A-Za-z0-9_ \\+\\-]+$"
    config <- read_yaml(info_file, fileEncoding = "UTF-8")
    if (is.null(config$author) || !(is.character(config$author)) || !(grepl(pattern, config$author))) {
        stop("Author is missing or has an incorrect format.")
    }
    organism <- tolower(config$organism)
    pattern <- "^(hg19|hg38|mm9|mm10|[a-z]{2}\\d{1,2})$"
    if (is.null(config$organism) || !(grepl(pattern, organism))) {
        stop("Organism is missing or has an incorrect format.")
    }
    vpchr <- as.character(config$VPchr)
    pattern <- "^(\\d{1,2}|X|Y)$"
    if (is.null(config$VPchr) || !(grepl(pattern, vpchr, ignore.case = TRUE))) {
        stop("VPchr is missing or has an incorrect format.")
    }
    if (is.null(config$VPpos) || !(is.numeric(config$VPpos))) {
        stop("VPpos is missing or has an incorrect format.")
    }
    check_reps <- function(reps) {
        if (!is.list(reps) && !is.vector(reps) || length(reps) < 2) {
            return(FALSE)
        }
        if (!all(sapply(reps, is.numeric))) {
            return(FALSE)
        }
        return(TRUE)
    }
    if (is.null(config$control)) {
        warning("control is missing.")
    }
    if (is.null(config$controlRep) || !(check_reps(config$controlRep))) {
        warning("controlRep is missing or has an incorrect format.")
    }
    if (!is.null(config$condition)) {
        if (is.null(config$conditionRep) || !(check_reps(config$conditionRep))) {
            stop("condition is provided, but conditionRep is missing or has an incorrect format.")
        }
    }
    if (is.null(config$readLength) || !(is.numeric(config$readLength))) {
        stop("readLength is missing or has an incorrect format.")
    }
    primer_f <- config$PrimerF
    primer_r <- config$PrimerR
    pattern <- "^[GATCatgc]+$"
    if (is.null(primer_f) || !(grepl(pattern, primer_f)) || is.null(primer_r) || !(grepl(pattern, primer_r))) {
        stop("One or both primer are missing or have an incorrect format.")
    }
    if (is.null(config$REEnz) || !(is.list(config$REEnz)) && !(is.vector(config$REEnz))) {
        stop("REEnz is missing or not a list.")
    }
    if (!all(sapply(config$REEnz, is.character))) {
        stop("REEnz must be a list where all entries are strings.")
    }
    if (!is.null(config$RESeq) && (is.list(config$RESeq) || is.vector(config$RESeq))) {
        is_valid_seq <- all(sapply(config$RESeq, function(seq) grepl("^[GATCatgc]+$", seq)))
        if (!(is_valid_seq)) {
            stop("RESeq contains invalid symbols")
        }
    } else {
        stop("RESeq is missing or not a list.")
    }
    return(TRUE)
}

