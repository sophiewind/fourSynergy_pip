library(BSgenome.Hsapiens.UCSC.hg19)


makeReSiteFile <- function(cutterSequence, chosenGenome) {

  chromosomeNames = seqnames(chosenGenome)
  totalSites = NULL
  siteChromosomes = NULL

  # human: take first 23 chromosomes (from chr1 to chrX)
  for (i in 1:23) {#length(chromosomeNames)) {
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



sites = makeReSiteFile("CATG", Hsapiens)
write.table(sites, file = "./scripts/fourSig/hg19_re_sites/b4_NlaIII_sites.txt", quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")


sites = makeReSiteFile("GATC", Hsapiens)
write.table(sites, file = "./scripts/fourSig/hg19_re_sites/b4_DpnII_sites.txt", quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")


sites = makeReSiteFile("AGATCT", Hsapiens)
write.table(sites, file = "./scripts/fourSig/hg19_re_sites/b6_BglII_sites.txt", quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")



sites = makeReSiteFile("GTAC", Hsapiens)
write.table(sites, file = "./scripts/fourSig/hg19_re_sites/b4_Csp6I_sites.txt", quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")
