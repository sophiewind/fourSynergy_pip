suppressMessages(library(yaml))


# Setup package ---------------------------------------------------------------
setwd('./scripts/fourSig/')
source('./fourSig.R')
setwd('../..')

#config <- read_yaml(file = './Datasets/Akincilar/info.yml')

# Receive args -----------------------------------------------------------------
tab <- commandArgs(trailingOnly = T)[1]  # TODO args
config <- read_yaml(file = commandArgs(trailingOnly = T)[2])  # TODO dynamic
OUT.PATH <- paste0("./results/",config$author, "/fourSig/")


# Setup directories ------------------------------------------------------------
# Create directory structure if it doesn't exist     
create_dir_if_not_exist <- function(dir_path) {
  if (!file.exists(dir_path)) {
    message('Builing dirs')
    dir.create(dir_path)
  }
}

create_dir_if_not_exist(paste0("./results/",config$author))
create_dir_if_not_exist(OUT.PATH)


for (n in c(1,3,5,11)){
  # Run analyses ---------------------------------------------------------------
  data <- fourSig(filename = tab,
                  chr = config$VPchr,
                  cis.only = T,
                  window.size = n,
                  iterations = 1000,
                  fdr = 0.01,
                  fdr.prob = 0.05)

  data.w <- readsToWindows(file = tab,
                           chr = config$VPchr,
                           window.size = n,
                           only.mappable = FALSE,
                           mask.start = config$VPpos-1000,  # TODO dyn
                           mask.end = config$VPpos+1000)

  # Plot results ---------------------------------------------------------------
  sample <- gsub('.tab', '', basename(tab))
  message(paste0(OUT.PATH, sample, '_fourSigit_', n, '.rds'))
  pdf(paste0(OUT.PATH, sample, '_plotit_', n, '.pdf'))
  try(plotReads(window.data = data.w,
                sig.data = data,
                x.start=config$VPpos-500000,
                x.end=config$VPpos+500000,  # TODO dym
                connect.points=TRUE,
                #mask.start = config$VPpos-5000,
                #mask.end = config$VPpos+5000,
                draw.reads.first=TRUE))
  dev.off()

  # Save rds -------------------------------------------------------------------

  saveRDS(data, paste0(OUT.PATH, sample, '_fourSigit_', n, '.rds'))
  saveRDS(data.w, paste0(OUT.PATH, sample, '_w_fourSigit_', n, '.rds'))


}
