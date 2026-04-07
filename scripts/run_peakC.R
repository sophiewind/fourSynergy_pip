library(devtools)
install_github("deWitLab/peakC")
library(peakC)
library(isotone)
library(caTools)
suppressMessages(library(yaml))

# Input from snakemake ---------------------------------------------------------
#wig.path <- dirname(commandArgs(trailingOnly = T)[1])
config <- read_yaml(file = commandArgs(trailingOnly = T)[1])  # TODO dynamic

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
# Process data with replicates
if (length(config$conditionRep) > 1){

  # Read data
  data.cond <- readMultipleWig(wigs.cond, vp.pos = viewpoint)
  if (!is.null(config$control)) data.ctrl <- readMultipleWig(wigs.ctrl, vp.pos = viewpoint)

  # Analyze data with different window sizes
  for (w in window.sizes){
    res.cond <- combined.analysis(data = data.cond, vp.pos = viewpoint, wSize = w)
    if (!is.null(config$control)) res.ctrl <- combined.analysis(data.ctrl, vp.pos = viewpoint, wSize = w)

    # Less stringent FDR
    res.cond.2 <- combined.analysis(data = data.cond, vp.pos = viewpoint, alphaFDR = 0.15, wSize = w)
    if (!is.null(config$control)) res.ctrl.2 <- combined.analysis(data.ctrl, vp.pos = viewpoint, alphaFDR = 0.15, wSize = w)

    # Visualize data and save .rds
    pdf(paste0(OUT.PATH,"peakC_", w, "_", config$condition, ".pdf"))
    plot_C(res.cond)
    dev.off()
    saveRDS(res.cond, paste0(OUT.PATH,"peakC_", w, "_", config$condition, ".rds"))

    if (!is.null(config$control)){
      pdf(paste0(OUT.PATH,"peakC_", w, "_", config$control, ".pdf"))
      plot_C(res.ctrl)
      dev.off()
      saveRDS(res.ctrl, paste0(OUT.PATH,"peakC_", w, "_", config$control, ".rds"))
    }

    # Less stringent FDR
    pdf(paste0(OUT.PATH,"peakCFDR_", w, " ", config$condition, ".pdf"))
    plot_C(res.cond.2)
    dev.off()
    saveRDS(res.cond.2, paste0(OUT.PATH,"peakCFDR_", w, " ", config$condition, ".rds"))

    if (!is.null(config$control)){
      pdf(paste0(OUT.PATH,"peakCFDR_", w, "_", config$control, ".pdf"))
      plot_C(res.ctrl.2)
      dev.off()
      saveRDS(res.ctrl.2, paste0(OUT.PATH,"peakCFDR_", w, "_", config$control, ".rds"))
    }
  }

# Process data without replicates  TODO rm
} else {
  # # Read data
  data.cond <- readqWig(wigs.cond, vp.pos = viewpoint, window = 700e3)
  if (!is.null(config$control)) data.ctrl <- readqWig(wigs.ctrl, vp.pos = viewpoint, window = 700e3)

  for (w in window.sizes){

    # Analyze data
    res.cond <- single.analysis(data = data.cond$data, vp.pos = viewpoint, qWd = 2.5, wSize = w)
    if (!is.null(config$control)) res.ctrl <- single.analysis(data.ctrl$data, vp.pos = viewpoint, qWd = 2.5, wSize = w)

    # Less stringent FDR
    res.cond.2 <- combined.analysis(data = data.cond, vp.pos = viewpoint, alphaFDR = 0.15, wSize = w)
    if (!is.null(config$control)) res.ctrl.2 <- combined.analysis(data.ctrl, vp.pos = viewpoint, alphaFDR = 0.15, wSize = w)

    # Visualize data and save .rds
    pdf(paste0(OUT.PATH,"peakC_", w, "_", config$condition, ".pdf"))
    plot_C(res.cond)
    dev.off()
    saveRDS(res.cond, paste0(OUT.PATH, config$condition,"peakC_", w, ".rds"))

    if (!is.null(config$control)){
      pdf(paste0(OUT.PATH,"peakC_", w, "_", config$control, ".pdf"))
      plot_C(res.ctrl)
      dev.off()
      saveRDS(res.ctrl, paste0(OUT.PATH,"peakC_", w, "_", config$control, "rds"))
    }

    # Less stringent FDR
    pdf(paste0(OUT.PATH,"peakCFDR_", w, "_", config$condition, ".pdf"))
    plot_C(res.cond.2)
    dev.off()
    saveRDS(res.cond.2, paste0(OUT.PATH,"peakCFDR_", w, "_", config$condition, ".rds"))

    if (!is.null(config$control)){
      pdf(paste0(OUT.PATH,"peakCFDR_", w, "_", config$control, ".pdf"))
      plot_C(res.ctrl.2)
      dev.off()
      saveRDS(res.ctrl.2, paste0(OUT.PATH,"peakCFDR_", w, "_", config$control, ".rds"))
    }
  }
}




