# Setup ------------------------------------------------------------------------
suppressMessages(library(r3Cseq))
#library(BSgenome.Mmusculus.UCSC.mm9.masked)  # TODO
suppressMessages(library(BSgenome.Mmusculus.UCSC.mm10.masked))  # dyn
suppressMessages(library(chromstaR))
suppressMessages(library(yaml))


#source('./FunctionsForBatchAnalysis_edited.R')


# Receive args -----------------------------------------------------------------
bam.path <- commandArgs(trailingOnly = T)[1]
config <- read_yaml(file = commandArgs(trailingOnly = T)[2])
OUT.PATH <- paste0("./results/",config$author, "/r3cseq/")

#config <- read_yaml(file = './Datasets/Fan/info.yaml')

# Setup directories ------------------------------------------------------------
# Create directory structure if it doesn't exist
create_dir_if_not_exist <- function(dir_path) {
  if (!file.exists(dir_path)) {
    message('Builing dirs')
    dir.create(dir_path)
  }
}

#create_dir_if_not_exist(paste0("./results/IMTB/"))
create_dir_if_not_exist(OUT.PATH)


# Read bams
# TODO pass this as param
bams <- list.files(bam.path, '*.bam$', full.names = F)
setwd(bam.path)
getwd()
if(!is.null(config$control)){
  message('\nControl involved\n')
  my3_obj <- new("r3CseqInBatch",
                 organismName=config$organism,  # TODO catch invalids
                 isControlInvolved=TRUE,  # TODO dyn
                 viewpoint_chromosome=paste0('chr', config$VPchr),
                 viewpoint_primer_forward=config$PrimerF,
                 viewpoint_primer_reverse=config$PrimerR,
                 restrictionEnzyme=config$REEnz[1],
                 bamFilesDirectory = '.',  # TODO möglich hier anzupassen?
                 BamExpFiles = bams[grep(config$condition, bams)],
                 BamContrFiles = bams[grep(config$control, bams)],
                 expBatchLabel = paste0(rep(config$condition, length(config$conditionRep)), '_', config$conditionRep),
                 contrBatchLabel = paste0(rep(config$control, length(config$controlRep)), '_', config$controlRep))
  
  # TODO: wants to write in read only file :( --> bad in docker 
  #setwd('../..')
  getwd()
  getBatchRawReads(my3_obj)
  
  
  # Get counts per fragment
  getBatchReadCountPerRestrictionFragment(my3_obj)
  
  
  # Normalize
  try(calculateBatchRPM(my3_obj))
  
  try(getBatchInteractions(my3_obj, method = 'intersection'))  # writes out interactions
  
  # visualize
  setwd('../../..')
  getwd()
  pdf(paste0(OUT.PATH, 'interactions_near_VP_r3Seq.pdf'))
  try(plotInteractionsNearViewpoint (my3_obj))
  dev.off()
  
  pdf(paste0(OUT.PATH, 'interactions_chr', config$VPchr,'_r3Seq.pdf'))
  try(plotInteractionsPerChromosome(my3_obj, paste0('chr', config$VPchr)))
  dev.off()
  
  pdf(paste0(OUT.PATH, 'interactions_overview_interactions_r3Seq.pdf'))
  try(plotOverviewInteractions(my3_obj))
  dev.off()
  #plotDomainogramNearViewpoint(my3_obj)
  
  
  try(saveRDS(my3_obj, paste0(OUT.PATH, 'r3Cseq_obj.rds')))
  
  
  
  
} else {
  message('\nNo control involved\n')
  message(getwd())
  tryCatch({
    for (bam in bams){
      message(bam)
      my3_obj <- new("r3Cseq", organismName=config$organism,
                     alignedReadsBamExpFile=bam,
                     isControlInvolved=FALSE, 
                     viewpoint_chromosome=paste0('chr', config$VPchr),
                     viewpoint_primer_forward=config$PrimerF,
                     viewpoint_primer_reverse=config$PrimerR,
                     expLabel= gsub('\"', "", deparse(bam, control = "all")),
                     restrictionEnzyme=config$REEnz[1])
      try(getRawReads(my3_obj))
      try(getReadCountPerRestrictionFragment(my3_obj))
      try(calculateRPM(my3_obj))
      try(getInteractions(my3_obj))  # writes out interactions
      
      setwd('../../..')
      getwd()
      pdf(paste0(OUT.PATH, gsub('.bam', '', (gsub('\"', "", deparse(bam, control = "all")))),
                 '_interactions_near_VP_r3Seq.pdf'))
      try(plotInteractionsNearViewpoint (my3_obj))
      dev.off()
      
      pdf(paste0(OUT.PATH, gsub('.bam', '', (gsub('\"', "", deparse(bam, control = "all")))), 'interactions_chr', config$VPchr,'_r3Seq.pdf'))
      try(plotInteractionsPerChromosome(my3_obj, paste0('chr', config$VPchr)))
      dev.off()
      
      pdf(paste0(OUT.PATH,gsub('.bam', '', (gsub('\"', "", deparse(bam, control = "all")))), 'interactions_overview_interactions_r3Seq.pdf'))
      try(plotOverviewInteractions(my3_obj))
      dev.off()
      #plotDomainogramNearViewpoint(my3_obj)
      
      message(paste0(OUT.PATH, gsub('.bam', '', basename(bam)), '_r3Cseq_obj.rds'))
      #try(
      saveRDS(my3_obj, paste0(OUT.PATH, gsub('.bam', '', basename(bam)), '_r3Cseq_obj.rds'))
      #)
     setwd(bam.path) 
    }
    
  }, error = function(e) {
      # Code to handle the error, for example:
      print("An error occurred but code execution continues")

})
}


