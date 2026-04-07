library(devtools)
install_github("deWitLab/peakC")  # TODO
library(peakC)
library(isotone)
library(caTools)
suppressMessages(library(yaml))



# Input from snakemake ---------------------------------------------------------
#wig.path <- dirname(commandArgs(trailingOnly = T)[1])
config <- read_yaml(file = commandArgs(trailingOnly = T)[1])

# Statics ----------------------------------------------------------------------
window.sizes <- c(1, 11, 21, 31, 51, 101)
OUT.PATH <- paste0("./results/",config$author, "/peakC/")

# Create directory structure if it doesn't exist
create_dir_if_not_exist <- function(dir_path) {
  if (!file.exists(dir_path)) {
    message('Builing dirs')
    dir.create(dir_path)
  }
}

create_dir_if_not_exist(paste0("./results/",config$author))
create_dir_if_not_exist(OUT.PATH)

# Define viewpoint
viewpoint <- config$VPpos  # Begin or end?
message(getwd())
#setwd('~')
#config <- read_yaml(file = './Datasets/IMTB/info.yml')
wigs.cond <- list.files(path=paste0("./results/", config$author, "/basic4cseq/wig"),
                 pattern=paste0('.*', config$condition, '.*'), full=T)
print(wigs.cond)
wigs.ctrl <- dir(path=paste0("./results/" ,config$author, "/basic4cseq/wig"), pattern =paste0('.*', config$control, '.*'), full=T)

# PeakC analysis ---------------------------------------------------------------
# Identifying interaction peaks
# Define the alphaFDR values to use
alphaFDR_values <- c(0.5, 0.25, 0.1, 0.05, 0.01)

# Read data
data.cond <- readMultipleWig(wigs.cond, vp.pos = viewpoint)
if (!is.null(config$control)) {
  data.ctrl <- readMultipleWig(wigs.ctrl, vp.pos = viewpoint)
}

# Analyze data with different window sizes and alphaFDR parameters
for (w in window.sizes) {
  for (alphaFDR in alphaFDR_values) {
    res.cond <- combined.analysis(data = data.cond, vp.pos = viewpoint, wSize = w, alphaFDR = alphaFDR)

    if (!is.null(config$control)) {
      res.ctrl <- combined.analysis(data.ctrl, vp.pos = viewpoint, wSize = w, alphaFDR = alphaFDR)
    }

    # Visualize data and save .rds for the main alphaFDR analysis
    pdf(paste0(OUT.PATH, "peakC_", w, "_", config$condition, "_alphaFDR_", alphaFDR, ".pdf"))
    plot_C(res.cond)
    dev.off()
    saveRDS(res.cond, paste0(OUT.PATH, "peakC_", w, "_", config$condition, "_alphaFDR_", alphaFDR, ".rds"))

    if (!is.null(config$control)) {
      pdf(paste0(OUT.PATH, "peakC_", w, "_", config$control, "_alphaFDR_", alphaFDR, ".pdf"))
      plot_C(res.ctrl)
      dev.off()
      saveRDS(res.ctrl, paste0(OUT.PATH, "peakC_", w, "_", config$control, "_alphaFDR_", alphaFDR, ".rds"))
    }
  }
}

