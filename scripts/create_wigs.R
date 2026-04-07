#!/usr/bin/env Rscript
# Setup ------------------------------------------------------------------------
message('\n####################################')
message('#     Loading libraries            #')
message('####################################\n')

suppressMessages(library(Basic4Cseq))
suppressMessages(library(yaml))
source('./scripts/basic4cseq_changed_wig.R')
source('./scripts/ensemble_functions.R')

# Statics ----------------------------------------------------------------------
config.yml <- commandArgs(trailingOnly = T)[2]
config <- read_yaml(file = config.yml)  # TODO dynamic

OUT.PATH <- paste0('./results/', config$author, '/basic4cseq/')


# Setup dirs -------------------------------------------------------------------
# Create directory structure if it doesn't exist
create_dir_if_not_exist <- function(dir_path) {
    if (!file.exists(dir_path)) {
        dir.create(dir_path)
    }
}

create_dir_if_not_exist(paste0("./results/",  config$author))
create_dir_if_not_exist(paste0("./results/vfl/"))
create_dir_if_not_exist(OUT.PATH)
create_dir_if_not_exist(paste0(OUT.PATH, "stats"))
create_dir_if_not_exist(paste0(OUT.PATH, "plots"))
create_dir_if_not_exist(paste0(OUT.PATH, "wig"))


# Debug ------------------------------------------------------------------------
message('\n####################################')
message('#     Starting Basic4CSeq          #')
message('####################################\n')
message(getwd())


# bam.path <- "../Einarbeitung_4C/example/bams/"
# config.yml <- "../rstudio/Datasets/IMTB/info.yml"



# Input from snakemake ---------------------------------------------------------
bam <- commandArgs(trailingOnly = T)[1]

#config <- read_yaml(file = './Datasets/IMTB/info.yml')
library(paste0('BSgenome.', ifelse(startsWith(config$organism, 'mm'), 'Mmusculus', 'Hsapiens'), '.UCSC.', config$organism, '.masked'), character.only = T)


# Input from config ------------------------------------------------------------
VP.chrom <- paste0('chr', config$VPchr)
ref <- config$organism  # TODO check and translate
re.1 <- config$RESeq[1]  # TODO check and translate
re.2 <- config$RESeq[2]
read.length <- determine_read_length_from_fastq(paste0('./Datasets/', config$author, '/',
                                                       gsub('_sorted.bam', '', basename(bam)),
                                                       '.fastq')) #ifelse(is.null(config$readLength), 34, config$readLength)


# Run Analysis -----------------------------------------------------------------
# TODO dynamic
points <- data.frame(chr = c(VP.chrom),
                     start = config$VPpos-1000,
                     end = config$VPpos+1000,
                     name = c('VP'),
                     colour = c('red'))


# Run Basic4Cseq   # TODO rm loop
# bam <- "/host/projects/swind/Einarbeitung_4C/example/bams/BMM_1_sorted_chr.bam"
print(bam)
ga <- readGAlignments(bam)
sample <- gsub('.bam', '', basename(bam))
sample <- gsub('_sorted', '', sample)

# Initialize object   TODO dynamic, für wig aber egal
data <- Data4Cseq(viewpointChromosome = VP.chrom, viewpointInterval = c(config$VPpos-1000, config$VPpos+1000),
                  readLength = read.length, pointsOfInterest = points, rawReads = ga)

# Get Fragments
# TODO dynamic
message(paste0('./results/vfl/', ref, '_', re.1, '_', re.2, '_', read.length, '_VFL.csv'))
if (!file.exists(paste0('./results/vfl/', ref, '_', re.1, '_', re.2, '_', read.length, '_VFL.csv'))){
    createVirtualFragmentLibrary_2(chosenGenome = get(paste0('BSgenome.', ifelse(startsWith(config$organism, 'mm'), 'Mmusculus', 'Hsapiens'), '.UCSC.', config$organism, '.masked')),
                                   firstCutter = re.1,
                                   secondCutter = re.2,
                                   readLength = read.length,
                                   libraryName = paste0('./results/vfl/', ref, '_', re.1, '_', re.2, '_', read.length, '_VFL.csv'))
}
fragmentLib = paste0('./results/vfl/', ref, '_', re.1, '_', re.2, '_',  read.length, '_VFL.csv')
print(fragmentLib)
rawFragments(data)<-readsToFragments_2(expData = data, fragmentLib = fragmentLib)
write.csv(rawFragments(data)[,c('chromosomeName', 'fragmentStart', 'fragmentEnd', 'fragEndReadsAverage')], file = paste0(OUT.PATH, sample, '_raw_frag.csv'),row.names = F, quote = F)
# Define near cis area
# TODO dynamic
if (nchar(config$RESeq[1]) == 6){
    nearCisFragments(data)<-chooseNearCisFragments(data,
                                                   regionCoordinates = c(config$VPpos-5000000, config$VPpos+5000000))
} else {
    nearCisFragments(data)<-chooseNearCisFragments(data,
                                                   regionCoordinates = c(config$VPpos-1000000, config$VPpos+1000000))
}

# Normalize Fragment data
nearCisFragments(data)<-normalizeFragmentData(data)
message(summary(nearCisFragments(data)))
# Get stats
getReadDistribution(data, outputName = paste0(OUT.PATH, 'stats/', sample, "_stats.txt"))
exportVisualizationFragmentData(data, fileName = paste0(OUT.PATH, sample, "_frags.csv"))

# Plot results
try(visualizeViewpoint(expData = data,  loessSpan = 0.05, mainColour = "blue", plotTitle = "", xAxisIntervalLength = 10000, maxY = 500))
try(drawHeatmap(expData = data, plotFileName = paste0(OUT.PATH, 'plots/', sample, "_heatmap.pdf"), xAxisIntervalLength = 10000))
try(visualizeViewpoint(data, plotFileName = paste0(OUT.PATH, 'plots/', sample, "_near_cis_max_1500.tiff"), loessSpan = 0.05, mainColour = "blue", plotTitle = "", xAxisIntervalLength = 10000, picDim = c(1800,1000), maxY = 1500))
try(drawHeatmap(data, plotFileName = paste0(OUT.PATH, 'plots/', sample, "_heatmap.tiff"), xAxisIntervalLength = 10000, picDim = c(1800,300)))

# Exctract VP chrom for wig file
data.chr <- data
data.chr@rawFragments <- data@rawFragments[data@rawFragments$chromosomeName == VP.chrom,]

# Create wig file
message('writing wig')
printWigFile2(data.chr, wigFileName = paste0("./results/", config$author, "/basic4cseq/wig/", sample, ".wig"))
message('end writing wig')
