library(rjson)
library(tibble)


prc.f1.collect <- data.frame(tool = c(), variable = c(), value = c(), dataset = c())
for (dataset in c('Geeven_sox', 'IMTB', 'Mermet', 'Jian',
                   # 'Akincilar', 
                  'Geeven_hba', 
                  #'Fan',
                  'Harris_c3_c4_VP+48','Harris_c3_c4_VP-20', 
                  'Harris_c12_c7_VP+48','Harris_c12_c7_VP-20',
                  'Mantis_Ldlrad4_TGFb',
                  'Hintermann_Hoxd1', 'Hintermann_Hoxd4', 'Hintermann_Hoxd9')){  #paste0('h4+4_DC', 2:4),  paste0('m4+4_DC', c(1,3,4)))
  f1 <- fromJSON(file = paste0('../results/', dataset, '/opo_single_44_prediction_f1_', dataset, '.json')) %>% as.data.frame() %>% t() %>% as.data.frame() 
  colnames(f1) <- 'weights f1 optimized'
  prc <- fromJSON(file = paste0('../results/', dataset, '/opo_single_44_prediction_prc_',  dataset, '.json')) %>% as.data.frame() %>% t() %>% as.data.frame()
  colnames(prc) <- 'weights prc optimized'
  f1.prc <- cbind(f1, prc)
  
  f1.prc$tool <- rownames(f1.prc)
  f1.prc.m <- reshape2::melt(f1.prc)
  f1.prc.m <-f1.prc.m# %>% 
    #mutate(value = ifelse(variable == 'weights f1 optimized', value, value *(-1)))
  f1.prc.m$dataset <- dataset
  prc.f1.collect <- rbind(prc.f1.collect, f1.prc.m)
}

#prc.f1.collect <-prc.f1.collect %>% mutate(origin = ifelse(startsWith(prefix =  'h', dataset), 'sim', 'public'))
#prc.f1.collect$dataset <- factor(prc.f1.collect$dataset,)  # c('Mermet', 'Geeven_sox', 'Geeven_hba', 'IMTB', paste0('h4+4_DC', 2:4),  paste0('m4+4_DC', c(1,3,4)))
ggplot(prc.f1.collect, aes(x = tool, y = value, fill = variable)) +
  geom_bar(stat = 'identity', position=position_dodge()) +
  facet_grid(dataset~.) +
  theme(axis.text.x = element_text(angle = -270))


prc.flat <- prc.f1.collect %>% 
  reshape2::dcast(., tool + dataset ~ variable, value.var = "value")


wilcox.test(prc.flat$`weights f1 optimized`, prc.flat$`weights prc optimized`, paired = TRUE)


# Compare ####
prf.collect <- data.frame(condition = c(), precision = c(), recall = c(),
                          f1 = c(), mcc = c(), auc_pr = c(),dataset = c())  #
for (dataset in c('Geeven_sox', 'IMTB', 'Mermet', 'Jian',
                  # 'Akincilar', 'Geeven_hba', 'Fan',
                  'Harris_c3_c4_VP+48','Harris_c3_c4_VP-20', 
                  'Harris_c12_c7_VP+48','Harris_c12_c7_VP-20',
                  'Mantis_Ldlrad4_TGFb',
                  'Hintermann_Hoxd1', 'Hintermann_Hoxd4', 'Hintermann_Hoxd9')){ #, 'Mermet', paste0('h4+4_DC', 2:4),  paste0('m4+4_DC', c(1,3,4)))){#, paste0('h4+4_DC', 2:4),  paste0('m4+4_DC', c(1,3,4)))){
 prf <- read.csv(paste0('../results/', dataset, '/precision_recall_f1_post.csv'))  %>% 
   mutate(method = ifelse(condition == 'known', 
                          'known', 
                          str_extract(condition, "(^[^.]+\\.[^.]+)|maj2..")))
 prf$dataset <- dataset
 
 prf.collect <- rbind(prf.collect, prf)
}


prf.collect.m <- melt(prf.collect)
prf.collect.m
prf.collect.m[is.na(prf.collect.m)] <- 0

#
prf.collect.m <- prf.collect.m[!grepl('single.*', prf.collect.m$condition), ]
prf.collect.m <- na.omit(prf.collect.m)

prf.collect.m <- prf.collect.m %>% 
  mutate(operation = sub("\\..*", "", condition)) %>%
  mutate(method = gsub('inter\\.', 'Intersection ', method)) %>% 
  mutate(method = gsub('rep\\.', '', method)) %>%
  mutate(method = gsub('\\+', ' and ', method)) %>%
  mutate(method = gsub('maj\\.', 'Majority vote ', method)) %>%
  mutate(method = gsub('union\\.', 'union ', method)) %>%
  mutate(operation = gsub('inter', 'intersection', operation)) %>%
  mutate(method = gsub('^union', 'Union ', method)) %>%
  mutate(operation = gsub('rep', 'tool', operation)) %>%
  #mutate(operation = gsub('\\+', ' and ', operation)) %>%
  mutate(operation = gsub('maj', 'majority vote', operation)) %>%
  mutate(method = gsub('r3c', 'r3Cseq', method)) %>%
  
  mutate(method = gsub('r4cker', '4C-ker', method)) %>%
  mutate(method = gsub('_', ' ', method)) %>%
  mutate(method = gsub('foursig', 'fourSig', method)) %>%
  mutate(method = gsub('fourcseq', 'FourCSeq', method)) %>%
  mutate(method = gsub('peakc', 'peakC', method)) %>%
  
  mutate(operation = gsub('double\\_', 'double ', operation)) %>%
  mutate(operation = factor(operation, levels = c('tool', 'union', 'intersection', 'majority vote'))) %>%
  #mutate(method = fct_reorder(method, as.numeric(operation), .desc = T)) %>% 
  mutate(variable = str_to_title(variable)) %>%
  mutate(variable = gsub('F1', 'F1 score', variable))
#

top <- as.numeric((prf.collect.m %>% group_by(condition, variable) %>% summarise(x = median(value)) %>% filter(variable == 'f1') %>% arrange(x) %>% tail(1))[, 'x'])

p <- ggplot(prf.collect.m, aes(y = method, x = value, fill = operation)) +  # $method, pattern = '(single|rep|bayes).*'),]  [prf.collect.m$dataset == 'Mermet',]
  geom_boxplot(alpha = 0.8) +  #linetype = freq  aes(fill = condition)
  geom_point() +
  facet_grid(~variable) +
 geom_vline(aes(xintercept =  top))

p

p + scale_fill_manual(values = levels.vec, guide = 'none') +
  scale_color_manual(values = c(all = 'black', missed = '#41333a'))

prf.collect.m %>% group_by(condition, variable) %>% summarise(x = mean(value)) %>% filter(variable == 'f1') %>% arrange(x) %>% head(1) %>% select(x)
  



# Test clustering
# todo center around vp
# Extract nearbait area
t1 <- read.csv('../results/Geeven_sox/basic4cseq/ESC_1_frags.csv', sep = '\t')
t2 <- read.csv('../results/Geeven_hba//basic4cseq/ESC_1_frags.csv', sep = '\t')

t3 <- read.csv('../results/Mermet/basic4cseq/KO_1_frags.csv', sep = '\t')
t4 <- read.csv('../results/Jian/basic4cseq/INS_A_beta_1_frags.csv', sep = '\t')


t5 <- read.csv('../results/Akincilar/basic4cseq/BLM14_1_frags.csv', sep = '\t')
t6 <- read.csv('../results/Court_Knck9_no///basic4cseq/D12_1_frags.csv', sep = '\t')

t7 <- read.csv('../results/Bian/basic4cseq/FL_1_frags.csv', sep = '\t')
t8 <- read.csv('../results/Vergult_VP5///basic4cseq/SHSY5Y_1_frags.csv', sep = '\t')

# Map on Fragments

# Cluster
counts.1 <- t1$reads
counts.2 <- t2$reads

counts.3 <- t3$reads
counts.4 <- t4$reads

counts.5 <- t5$reads
counts.6 <- t6$reads

# TODO normalize read counts


df <- (data.frame(gs = counts.1[1:1000], gh = counts.2[1:1000], 
                 mem = counts.3[1:1000], jian = counts.4[1:1000],
                 bian = counts.5[1:1000], vg5 = counts.6[1:1000]))

ComplexHeatmap::Heatmap(df, cluster_columns = F)

# Remove 0 var columns
df <- t(data.frame(gs = counts.1[1:1000], gh = counts.2[1:1000], 
                  mem = counts.3[1:1000], jian = counts.4[1:1000],
                  bian = counts.5[1:1000], vg5 = counts.6[1:1000]))
variances <- apply(df, 2, var)
df.filtered <- df[, variances > 0]

pca_result <- prcomp(df)
plot(pca_result$x[,1:2], pch = 19, xlab = "PC1", ylab = "PC2")


dist_matrix <- dist(df)
# Hierarchisches Clustering durchführen
hclust_result <- hclust(dist_matrix, method = "complete") 
plot(hclust_result)


# pCA
log_counts <- log2(df.filtered + 1)
dataset_labels <- factor(c("Dataset1", "Dataset2", "Dataset3", "Dataset4", "Dataset5", "Dataset6"))
# Durchführung der PCA
pca_result <- prcomp(log_counts, center = TRUE, scale. = TRUE)

# Ergebnisse anzeigen
summary(pca_result)  # Zeigt die Proportion der erklärten Varianz
print(pca_result$rotation)  # Zeigt die Variablenladungen

pca_scores <- as.data.frame(pca_result$x)
pca_scores$Dataset <- dataset_labels
ggplot(pca_scores, aes(x = PC1, y = PC2, color = Dataset, label = Dataset)) +
  geom_point(size=4) +
  # geom_text(vjust=-1, hjust=0.5, size=3) +
  labs(title="PCA of Datasets", x="Principal Component 1", y="Principal Component 2") +
  theme_minimal() 
