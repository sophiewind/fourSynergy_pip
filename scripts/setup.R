library(devtools)
install_github("deWitLab/peakC")
BiocManager::install("r3Cseq")  # update none
BiocManager::install("BSgenome.Mmusculus.UCSC.mm9.masked")  # update none
BiocManager::install("BSgenome.Mmusculus.UCSC.mm10.masked")
BiocManager::install("BSgenome.Mmusculus.UCSC.mm39.masked")
BiocManager::install("BSgenome.Mmusculus.UCSC.mm39")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg19.masked")
BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene")
BiocManager::install('DESeq2')
BiocManager::install("chromstaR")
BiocManager::install("Basic4Cseq")

install.packages('isotone')  # peakC
install.packages('caTools')  # peakC
install.packages('UpSetR')

install.packages('tidyverse')

install.packages('cowplot')
install.packages('patchwork')


install.packages('gridExtra')

install.packages('valr')

BiocManager::install("ggbio")
# MArkdown

install.packages('ggbreak')

# install.packages('renv')
# renv::init()

# conda install bioconda::deeptools

install.packages('svglite')

BiocManager::install("karyoploteR")
