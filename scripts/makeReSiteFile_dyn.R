library(yaml, quietly = T)

config.yml <- paste0(commandArgs(trailingOnly = T)[1])


config <- read_yaml(file = config.yml)
library(paste0('BSgenome.', ifelse(startsWith(config$organism, 'mm'), 'Mmusculus', 'Hsapiens'), '.UCSC.', config$organism, '.masked'), character.only = T, quietly = T)


makeReSiteFile <- function(cutterSequence, chosenGenome) {
    chromosomeNames = seqnames(chosenGenome)
    totalSites = NULL
    siteChromosomes = NULL

    # human: take first 23 chromosomes (from chr1 to chrX)
    for (i in 1:ifelse(config$organism == 'hg19', 23, 21)){#length(chromosomeNames)) {
        print(i)
        chromosomeToSplit = chosenGenome[[i]]
        # subtract 1 from position to match original fourSig file structure
        currentSites = digestChromosome(cutterSequence, chromosomeToSplit) - 1
        totalSites = c(totalSites, currentSites)
        chrName = gsub("chr", "", seqnames(chosenGenome)[i])
        siteChromosomes = c(siteChromosomes, rep(chrName, length(currentSites)))
    }

    restrictionSites = data.frame("chr" = siteChromosomes, "pos" = totalSites)
    return(restrictionSites)
}

digestChromosome <- function(cutterSequence, chromosomeToSplit) {
    dnaSequence = toString(as(unmasked(chromosomeToSplit), "Views"))
    fragmentSequences = gregexpr(toupper(cutterSequence), dnaSequence)[[1]]
    return(fragmentSequences)
}


# Get right genome
if (config$organism == 'hg19') {
    genome <- Hsapiens
} else if (config$organism == 'mm10') {
    genome <- Mmusculus
} else {
    stop("Unsupported organism specified in config")
}

if (!dir.exists(paste0("./scripts/fourSig/", config$organism, "_re_sites/"))){
    dir.create(paste0("./scripts/fourSig/", config$organism, "_re_sites/"))
}

# Now you can use the 'genome' object
print(genome@pkgname)

if (!file.exists(paste0("./scripts/fourSig/", config$organism, "_re_sites/b", nchar(config$RESeq[1]), "_", config$REEnz[1], "_sites.txt"))){
    message(paste0("creating ./scripts/fourSig/", config$organism, "_re_sites/b", nchar(config$RESeq[1]), "_", config$REEnz[1], "_sites.txt\n"))
    sites = makeReSiteFile(config$RESeq[1], genome)
    write.table(sites, file = paste0("./scripts/fourSig/", config$organism, "_re_sites/b", nchar(config$RESeq[1]), "_", config$REEnz[1], "_sites.txt"), quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")
}

if (!file.exists(paste0("./scripts/fourSig/", config$organism, "_re_sites/b", nchar(config$RESeq[2]), "_", config$REEnz[2], "_sites.txt"))){
    message(paste0("creating ./scripts/fourSig/", config$organism, "_re_sites/b", nchar(config$RESeq[2]), "_", config$REEnz[2], "_sites.txt\n"))
    sites = makeReSiteFile(config$RESeq[2], genome)
    write.table(sites, file = paste0("./scripts/fourSig/", config$organism, "_re_sites/b", nchar(config$RESeq[2]), "_", config$REEnz[2], "_sites.txt"), quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")
}




