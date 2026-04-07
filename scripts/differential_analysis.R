library(DESeq2)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
config <- read_yaml('./Datasets/Geeven_hba////info.yaml')
setwd("~/work/")
res.path <- paste0('./results/', config$author, '/')

# TODO sind diese counts normalized, da keine ints? ?Craolin TODO
# Get read counts per fragment from Basic
basic.frag <- read.delim(paste0(res.path, '/nearbait_area.bed'),  sep = '\t')  # TODOdyn

for (b in list.files(paste0(res.path, 'basic4cseq/'), pattern = 'raw_frag.csv', full.names = T)[-1]){
  print(b)
  sample <- gsub('_raw_frag.csv', '', basename(b))
  basic.frag.tmp <- read.csv(b)
  colnames(basic.frag.tmp)[4] <- paste0('reads_', sample)
  basic.frag <- merge(basic.frag, basic.frag.tmp, all = T)
}

ctrl <- basic.frag[grep(pattern = config$condition, colnames(basic.frag))]

# Mean counts between replicates (simes, love obsidian)
# basic.frag <- basic.frag %>%
#   mutate(!!sym(config$condition) := rowMeans(select(., contains(paste0('reads_', config$condition))), na.rm = TRUE)) %>%
#   mutate(!!sym(config$control) := rowMeans(select(., contains(paste0('reads_', config$control))), na.rm = TRUE)) 

# Remove raw data
# Count data
library(tibble)
#basic.frag <- basic.frag %>% select(-starts_with('reads_'))
colnames(basic.frag) <- gsub('reads_', '', colnames(basic.frag)) 
colnames(basic.frag)[1:3] <- c('seqnames', 'start', 'end')


basic.frag.gr <- makeGRangesFromDataFrame(basic.frag, keep.extra.columns = T)

# TODO normalize counts
# NO: unnormalized
# bio reps should not be collapse with collapseRep fct be DESeq2


# prepare colData
coldata <- data.frame('condition' = rep(c('condition', 'control'), each = length(config$conditionRep)))  # order should be definded bey mutate(?) 
rownames(coldata) <- names(mcols(basic.frag.gr))


# Nearbait ---------------------------------------------------------------------
window <- ifelse(nchar(config$RESeq)[1] == 4, 1000000, 5000000)
basic.frag.gr.nb <- subsetByOverlaps(basic.frag.gr, GRanges(paste0('chr', config$VPchr), IRanges(config$VPpos-window, config$VPpos+window)))  # TODO dyn
counts.nb <- basic.frag.gr.nb %>% 
  as.data.frame() %>% 
  mutate(pos = paste(seqnames(basic.frag.gr.nb), start(basic.frag.gr.nb), end(basic.frag.gr.nb), sep = ':')) %>%
  column_to_rownames(var = 'pos') %>%
  dplyr::select(-seqnames, -start, -end, -width, -strand)

dds <- DESeqDataSetFromMatrix(countData =counts.nb,  # TODO change, this cant be valid
                              colData = coldata,
                              design = ~ condition)
resultsNames(dds)
# Shrink?
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds, fitType = "local")
dds <- nbinomWaldTest(dds)
res <- results(dds, contrast = c("condition", 'condition', 'control'))  # correct?
norm_counts <- counts(dds, normalized = TRUE)
norm_counts_log <- log(norm_counts+1,10)

plotMA(res, ylim=c(-2,2))
plotCounts(dds, gene=which.min(res$padj), intgroup="condition")

#res.intervall
library(tidyr)
res.intervall.diff <- res %>% as.data.frame() %>% 
  rownames_to_column('pos') %>%
  separate(pos, c('seqnames', 'start', 'end')) %>% 
  filter(padj < 0.05)

res.intervall.diff.gr <- makeGRangesFromDataFrame(res.intervall.diff, keep.extra.columns = T)

library(karyoploteR)
pp <- getDefaultPlotParams(plot.type = 2)
pp$ideogramheight <- 10
pp$data2height <- 75
kp <-plotKaryotype(genome = config$organism, chr = paste0('chr',config$VPchr), zoom=toGRanges(paste0('chr',config$VPchr,':', min(start(res.intervall.diff.gr))-10000, ':',max(start(res.intervall.diff.gr))+10000)), plot.type = 2, plot.params = pp)

#kp <-plotKaryotype(genome = config$organism, chr = paste0('chr',config$VPchr), zoom=toGRanges(paste0('chr',config$VPchr,':',config$VPpos-10000 , ':',config$VPpos+10000)), plot.type = 2)

# genes.data <- makeGenesDataFromTxDb(TxDb.Hsapiens.UCSC.hg19.knownGene,  # TODO
#                                     karyoplot=kp,
#                                     plot.transcripts = TRUE, 
#                                     plot.transcripts.structure = TRUE)

# TODO collapse rep bigwig
#kp <-plotKaryotype(genome = config$organism, chr = paste0('chr',config$VPchr), zoom=toGRanges(paste0('chr',config$VPchr,':', min(start(res.intervall.diff.gr))-10000, ':',max(start(res.intervall.diff.gr))+10000)), plot.type = 2)
#kpPlotRegions(kp, res.intervall.diff.gr, r0 = 0.25, r1 = 0.75, col = 'lightgrey')
kpPoints(kp,chr = paste0('chr',config$VPchr), x = start(res.intervall.diff.gr), y = 0, pch = 17, cex = 2, col = 'black')
kpPlotBAMCoverage(kp, data=paste0('./results/', config$author, '/alignment/', config$condition,'_1_sorted.bam'), col="#1E90FF", ymax=1000, max.valid.region.size = 2000001, r0 = 0.5, r1 = 0)
kpPlotBAMCoverage(kp, data=paste0('./results/', config$author, '/alignment/', config$control,'_1_sorted.bam'), col="red", ymax=1000, max.valid.region.size = 2000001, r0=0.5,r1 = 1)
#kpPlotGenes(kp, data=genes.data)

#kpPoints(kp,chr = paste0('chr',config$VPchr), x = start(res.intervall.diff.gr), y = 0, pch = 17, cex = 2, col = 'black')
kpBars(kp, chr = paste0('chr',config$VPchr), x0=start(res.intervall.diff.gr), x1 = end(res.intervall.diff.gr), 
       y0 = 0, y1=log10(res.intervall.diff$padj)*-1, data.panel = 2,
       ymax = max(log10(res.intervall.diff$padj)*-1))
kpAddLabels(kp, 'coverage', label.margin = 0.035)
kpAddLabels(kp, '-log10(padj)', label.margin = 0.035, data.panel = 2)
#kpAxis(kp, data.panel=1, r0 = 0)
kpAxis(kp, data.panel=2, ymin = 0, ymax = max(log10(res.intervall.diff$padj)))
kpRect()
# ggpot
res.intervall.diff.gr

ggplot()

# volcano

# Cis --------------------------------------------------------------------------