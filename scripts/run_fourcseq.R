# setwd('/home/rstudio/data/')
#BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
library(FourCSeq)

suppressMessages(library(yaml))

#config <- read_yaml(file = '/home/rstudio/work//Datasets/Fan/info.yaml')

config <- read_yaml(file = commandArgs(trailingOnly = T)[2])


if (config$organism == 'mm10'){
  suppressMessages(library(TxDb.Mmusculus.UCSC.mm10.knownGene))
} else if (config$organism == 'hg19'){
  suppressMessages(library(TxDb.Hsapiens.UCSC.hg19.knownGene))
} else {
  message('no valid org (fourcseq)')
}

ref <- commandArgs(trailingOnly = T)[1]
#ref <- "/host/genomes/Mus_musculus.GRCm38.p6/bwa.0.7.17/GRCm38.p6.genome_clean.fa" 
#ref <- "/host/genomes/Homo_sapiens.GRCh37.67/bwa-0.7.17/chr/Homo_sapiens.GRCh37.67.dna.chromosome.all.chr.fasta"
writeTo <- OUT.PATH <- paste0("./results/",config$author, "/fourcseq/")


## ###########
## FourCSeq
## ###########

print("FourCSeq")

bamFilePath = paste0("./results/", config$author, 
                     "/alignment")

metaData <- list(projectPath = writeTo, fragmentDir = "re_fragments", 
                 referenceGenomeFile = ref, 
                 reSequence1 = config$RESeq[1],
                 reSequence2 = config$RESeq[2], 
                 primerFile = "", bamFilePath = bamFilePath)

# colData <- DataFrame(viewpoint = "Test", 
#                      condition = factor(c(rep(config$condition, length(config$condition)), rep(config$control, length(config$control))), 
#                                         levels = c(config$condition, config$control)),
#                      replicate = c(1:length(labelCondition),1:length(labelControl)), bamFile = gsub("_chr", "_sorted", c(bamsCondition, bamsControl)), sequencingPrimer = "first")
# colLevels = levels(colData$condition)
## TODO more elegant
if (is.null(config$control)) {
  bams <- c(paste0(config$condition,
                   "_",
                   config$conditionRep,
                   "_sorted.bam"))
} else {
  bams <- c(paste0(config$condition,
                   "_",
                   config$conditionRep,
                   "_sorted.bam"),
            paste0(config$control,
                   "_",
                   config$controlRep,
                   "_sorted.bam"))
}

colData <- DataFrame(viewpoint = "VP",
                     condition = factor(c(rep(config$condition, 
                                              length(config$conditionRep)),  # Achtung order!!! TODO test
                                          rep(config$control, 
                                              length(config$controlRep)))),
                     replicate = c(config$conditionRep, config$controlRep),
                     
                     # Check if control is involved
                     bamFile = bams,
                     sequencingPrimer="first")
colData

fc1 <- FourC(colData, metaData)

fc1 <- addFragments(fc1)

# manipulate fc1 via rowRanges here to insert valid near-cis area around the viewpoint for test run (GRanges object)!

colData(fc1)$chr = paste0('chr', config$VPchr)
colData(fc1)$start = config$VPpos
colData(fc1)$end = config$VPpos + 1000

fc1 <- countFragmentOverlaps(fc1, trim=4, minMapq=30)
fc1 <- combineFragEnds(fc1)
fc1 <- smoothCounts(fc1)

writeTrackFiles(fc1)

#, type = 'iterate'
try(fcf1 <- getZScores(fc1, minCount = 20, minDist = 1000)     )

try(fcf1 <- addPeaks(fcf1, zScoreThresh=1, fdrThresh=0.05))
#if (exists('fcf1'))
tryCatch({
 
  head(SummarizedExperiment::assay(fcf1, 'peaks'))
  table(assay(fcf1, 'peaks') )
  
  for (i in seq_along(colnames(fcf1))){
    sample.tmp <- gsub('VP_', '', colnames(fcf1)[i])
    
    # Extract peak range
    peaks.gr <- rowRanges(fcf1)[assay(fcf1, 'peaks')[,i],0]
    
    if (length(peaks.gr) > 0){
      peaks.gr$tool <- paste0('single.fourcseq.', sample.tmp)
    } else {
      peaks.gr <- makeGRangesFromDataFrame(data.frame(seqnames='chr1',start=0,end=0, nReads = -1, tool = paste0('single.fourcseq.', sample.tmp)), keep.extra.columns = T)
    }
    message(paste0('writing to: ',writeTo, 'sia_fourcseq_', sample.tmp, '.rds'))
    saveRDS(peaks.gr, paste0(writeTo, sample.tmp, '_fourcseq.rds'))
    # df <- df %>% filter(peak == 1)
    # try(df$seqnames <- paste0('chr', config$VPchr))  # TODO richtig?? nur cis)
    # try(df$tool <- paste0('single.fourcseq.', sample.tmp))
    # try(list.sia.all[[paste0('single.fourcseq.', sample)]] <- makeGRangesFromDataFrame(df, keep.extra.columns = T))  # TODO empty gr
  }
  
  
  if (!is.null(config$control)){
    fcf1 <- getDifferences(fcf1, config$condition)
    res <- getAllResults(fcf1)
    saveRDS(res, paste0(writeTo, 'res.rds'))
  }

  
  saveRDS(fcf1, paste0(OUT.PATH, 'fcf.rds'))
  saveRDS(list(
    assays = fcf1@assays,
    rowRanges = fcf1@rowRanges,
    colData = fcf1@colData), 
    paste0(OUT.PATH, 'fcf_list.rds'))
  
  
  png(paste0(writeTo, 'PlotZScore.png'))
  txdb <- if (config$organism =='mm10') {
    TxDb.Mmusculus.UCSC.mm10.knownGene
  } else {
    TxDb.Hsapiens.UCSC.hg19.knownGene
  }
  
  plotZScores(fcf1, txdb = txdb)
  dev.off()
}, error = function(e) {
  message("An error occurred:")
  print(e)
  # You can also write the error to a log file or perform other actions
}
)


