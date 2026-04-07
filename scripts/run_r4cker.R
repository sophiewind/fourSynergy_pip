# Setup ------------------------------------------------------------------------
#library(devtools)
library(R.4Cker)
library(ggplot2)
library('dplyr')
library(yaml)
#install_github("rr1859/R.4Cker")  # conda skeletom


# Receive args -----------------------------------------------------------------
# config <- read_yaml(file = './Datasets/geeven_sox/info.yaml')
config <- read_yaml(file = commandArgs(trailingOnly = T)[1])

# Setup directories ------------------------------------------------------------
OUT.PATH <- paste0('./results/', config$author, '/r4cker/')

# Create directory structure if it doesn't exist
create_dir_if_not_exist <- function(dir_path) {
  if (!file.exists(dir_path)) {
    message('Builing dirs')
    dir.create(dir_path)
  }
}

create_dir_if_not_exist(paste0("./results/",config$author))
create_dir_if_not_exist(OUT.PATH)


# Run oligoMatch ---------------------------------------------------------------
# Write RE sequence fasta if not already existing
if (!file.exists(paste0('./scripts/R.4Cker/', tolower(config$REEnz[1]), '.fa'))){
  fileConn<-file(paste0('./scripts/R.4Cker/', tolower(config$REEnz[1]), '.fa'))
  writeLines(text = config$RESeq[1])
  close(fileConn)
}

message((paste0('oligoMatch ./scripts/R.4Cker/', tolower(config$REEnz[1]), '.fa ',
                      config$ref, ' ',
                      OUT.PATH,tolower(config$REEnz[1]),'_',config$organism,'_restriction_side_oligomatch.bed')))
system(paste0('oligoMatch ./scripts/R.4Cker/', tolower(config$REEnz[1]), '.fa ',
              config$ref, ' ',
              OUT.PATH,tolower(config$REEnz[1]),'_',config$organism,'_restriction_side_oligomatch.bed'))
system(paste0('oligoMatch ./scripts/R.4Cker/', tolower(config$REEnz[1]), '.fa ',
              config$ref, ' ',
              OUT.PATH,tolower(config$REEnz[2]),'_',config$organism,'_restriction_side_oligomatch.bed'))


# R4cker analyses --------------------------------------------------------------
enz_file=read.table(paste0(OUT.PATH, tolower(config$REEnz[1]),'_',
                           config$organism,'_restriction_side_oligomatch.bed'),
                    stringsAsFactors = FALSE)

# Define number of replicates
if(is.null(config$control)) {
  rep <- max(config$conditionRep)
} else {
  reps <- c(max(config$conditionRep), max(config$controlRep))
}

# Create r4cker file
my_obj = createR4CkerObjectFromFiles(files = list.files(paste0(gsub('r4cker/', '', OUT.PATH), 'alignment/'), '*sorted_rm_self_und.bedGraph$', full.names = T),
                                     bait_chr=paste0('chr', config$VPchr),
                                     bait_coord= config$VPpos,
                                     bait_name = config$author,
                                     primary_enz = config$RESeq[1],
                                     samples = c(paste(config$condition, (config$conditionRep),sep =  '_'),
                                                 paste(config$control, (config$controlRep),sep =  '_')),  #TODO code softer
                                     conditions = c(config$condition, config$control),
                                     replicates = reps,
                                     species = ifelse(substr(config$organism, 1, 2) == "mm", "mm", "hs"),  # TODO eleganter
                                     output_dir = paste0(OUT.PATH, 'out'),
                                     enz_file=enz_file)


#' Output:
#' - QC statistics
#' - normalized counts for all together and indiv samples
#' - highinter for replicates together
#' - ??C high inter entscheident
#' - BMM (hi. low. non): 408; 81; 433 -> counts ca. 3.8k
#' ??Warum nur an einer Stelle counts
#' ??was ist adaptive Window


# ??VP nicht ausreichend gecovert?
nb_results=nearBaitAnalysis(my_obj,k=5)
ggplot(nb_results$norm_counts_avg, aes(x=Coord, y=Count, colour=Condition))+
  theme_bw()+
  geom_line()+xlab(paste("Chromosome coordinates (", my_obj@bait_chr, ")", sep =""))+
  ylab("Normalized counts")+
  ggtitle(paste("Near bait analysis (", my_obj@bait_name, " bait)", sep = ""))
ggsave(paste0(OUT.PATH, 'out/nearBait.pdf'))
# TODO rausschreiben, was Carolin auch rauschreibt


# nb_results_3=nearBaitAnalysis(my_obj,k=3)
# ggplot(nb_results_3$norm_counts_avg, aes(x=Coord, y=Count, colour=Condition))+
#   theme_bw()+
#   geom_line()+xlab(paste("Chromosome coordinates (", my_obj@bait_chr, ")", sep =""))+
#   ylab("Normalized counts")+
#   ggtitle(paste("Near bait analysis (", my_obj@bait_name, " bait)", sep = ""))

cis_results=cisAnalysis(my_obj,k=10)
ggplot(cis_results$norm_counts_avg, aes(x=Coord, y=Count, colour=Condition))+
  theme_bw()+geom_line()+xlab(paste("Chromosome coordinates (", my_obj@bait_chr, ")", sep =""))+
  ylab("Normalized counts")+
  ggtitle(paste("cis (", my_obj@bait_name, " bait)", sep = ""))
ggsave(paste0(OUT.PATH, 'out/cis.pdf'))


try(transAnalysis(my_obj,k=20))

# Differential Interactions
if (!is.null(config$control)){
  res_df = differentialAnalysis(obj=my_obj,
                                norm_counts_avg=nb_results$norm_counts_avg,
                                windows=nb_results$window_counts,
                                conditions=c(config$condition, config$control),
                                region="nearbait",
                                coordinates=NULL,
                                pval=0.05)

  write.csv(res_df, paste0(OUT.PATH, 'differentialRegions.csv'))

  res_df %>% mutate(negLog10 = -log10(padj)) %>%
    ggplot(., aes(x = log2FoldChange, y = negLog10)) +
    geom_point()
  ggsave(paste0(OUT.PATH, 'out/differential.pdf'))
}


