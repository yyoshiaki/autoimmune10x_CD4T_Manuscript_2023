library(tidyverse)
library(data.table)
library(ggbeeswarm)
library(pheatmap)
library(viridis)
library(RColorBrewer)
library(extrafont)
library(patchwork)

loadfonts()

setwd('/home/rstudio/autoimmune_10x')
prefix.output <- "./output/metaanalysis/230816_CD4T"
project.name.all <- c("GSE157278", "GSE135779", "GSE132338", "E-MTAB-10026",
                      "GSE168732", "GSE163314", "GSE194315", "CeD", "577860",
                      "GSE144744", "PRJNA605083", "ours", "3993994", "HRA000155",
                      "GSE174188", "GSE138266", "GSE149689", "GSE150728", 
                      "GSE185857", "GSE158055", "GSE133028", "GSE181279", "1M_PBMC_T1D_Parse", 
                      "GSE198616", "GSE125527")
# blacklist_disease <- c("Control_ParkinsonDisease", "Control_BirdshotChorioretinitis")
blacklist_disease <- c()
f_blacklist <- './data/230122_blacklist_samples.csv'

df.blacklist <- read_csv(f_blacklist)
blacklist_sample <- paste(df.blacklist$project, df.blacklist$sample, sep = ":")

### count
list.df <- list()

for (project.name in project.name.all) {
  df.metadata <- read_csv(paste0('data/', project.name,'/', project.name,'_metadata.csv'))
  df.output <- read_csv(paste0('output/', project.name, '/',  project.name, '_queryL2_Reference_Mapping.csv'))
  
  df.output <- df.output %>% 
    left_join(df.metadata, by = "sample") %>%
    group_by(sample) %>%
    mutate(freq = n / sum(n), project = project.name,
           sample = paste(project.name, sample, sep = ":"))
  
  list.df[[project.name]] <- df.output
}

df.output <- rbindlist(list.df, fill = TRUE) %>% 
  as_tibble() %>%
  mutate(sex = str_to_lower(sex), age = age / 25)

df.output <- df.output %>%
  filter(!disease %in% blacklist_disease) %>%
  filter(!sample %in% blacklist_sample)

df.output$disease <- str_replace(df.output$disease, "COVID-19", "Covid")
# 230112 CD&SpA to CD and droped SpA
df.output$disease <- str_replace(df.output$disease, "CD&SpA", "CD")
df.output <- df.output %>% filter(disease != "SpA") %>%
  mutate(disease = relevel(factor(disease), ref = "HC"))
write_csv(df.output, paste0(prefix.output, '_tidy_queryL2_output.csv'))

n_total <- df.output %>% group_by(project) %>% summarise(n_total = sum(n))
n_hc <- df.output %>% group_by(project) %>% filter(disease == "HC") %>% summarise(n_HC = n_distinct(sample))
n_dis <- df.output %>% group_by(project) %>% filter(disease != "HC") %>% summarise(n_dis = n_distinct(sample))
df.stats <- n_total %>% left_join(n_hc, by = "project") %>% left_join(n_dis, by = "project")
write_csv(df.stats, paste0(prefix.output, '_stats.csv'))

p <- ggplot(df.output, aes(x=clusterL2, y=freq, color=disease), scale = "width") + 
  geom_quasirandom(dodge.width=1, size=0.5) +
  theme(axis.text.x = element_text(angle = 90))
p
ggsave(paste0(prefix.output, '_queryL2_ggbeeswarm.pdf'), width = 10, height = 4)

df.ntotal <- df.output %>% group_by(sample) %>% summarise(n_total = sum(n))
cols <- c("sample", "disease", "disease_duration", "condition", "age", "race", 
"sex", "immunosuppressant", "immunosuppressant_duration", "project")
df.metadata <- df.output[,cols] %>% distinct()


df.output.wide <- df.output %>% 
  pivot_wider(values_from = n, names_from = clusterL2,
              id_cols = c(sample, disease, disease_duration, condition, age, race, sex,
                          immunosuppressant, immunosuppressant_duration, project), values_fill = 0) %>%
  left_join(df.ntotal, by='sample') %>%
  drop_na(age)

df.output.wide <- as.data.frame(df.output.wide)
df.output.wide <- within(df.output.wide, disease <- relevel(factor(disease), ref = "HC"))
df.output.wide <- within(df.output.wide, sex <- relevel(factor(sex), ref = "male"))

cells <- c("Tnaive", "TnaiveAct", "TnaiveMX1", "TnaiveSOX4","TcmTh0","TcmTh0Act", 
           "TcmTfh", "TcmPHLDA3", "TcmTh17", "TcmTh2",
           "TemTh1pre", "TemTh1", "TemTh117", "TemTph", "TemraTh1",
           "TregNaive", "TregAct", "TregEff")

list.res <- list()
for (cell in cells) {
  res <- glm(as.formula(paste0("cbind(", cell, ", n_total-", cell,") ~ disease + age + sex + project")),
             data = df.output.wide, family = binomial)
  # res <- glm(as.formula(paste0("cbind(", cell, ", n_total-", cell,") ~ disease + age + sex")),
  #            data = df.output.wide, family = binomial)
  list.res[[cell]] <- as_tibble(summary(res)$coefficients, rownames = "var") %>%
    mutate(cell = cell)
}

df.res <- rbindlist(list.res, fill = TRUE) %>% 
  as_tibble() %>%
  mutate(padj = p.adjust(`Pr(>|z|)`))

df.plot <- df.res %>%
  filter(str_detect(var, "^disease|^age|^sex")) %>%
  mutate(cell = factor(cell, levels = cells))

lim.est <- 1.8
lim.nest <- -1.8
lim.lpa <- 10^-100
lim.lpa.lower <- 0.05

df.plot <- df.plot %>% 
  mutate(Estimate=if_else(Estimate>lim.est, lim.est , Estimate)) %>% 
  mutate(Estimate=if_else(Estimate< lim.nest, -lim.est, Estimate)) %>%
  mutate(padj=if_else(padj<lim.lpa, lim.lpa, padj)) %>%
  mutate(var=fct_relevel(var, "sexfemale", "age"))

order.row <- df.plot %>% group_by(var) %>% filter(`z value` == max(`z value`)) %>% arrange(cell) %>% select(var)
order.row <- as.character(order.row$var)
order.row <- order.row[!order.row %in% c("age", "sexfemale")]
order.row <- append(order.row, c('age', 'sexfemale'))
# df.plot$var <- factor(df.plot$var, levels = rev(order.row))
df.plot$var <- str_replace(df.plot$var, 'disease', 'dis.')
df.plot$var <- factor(df.plot$var, levels = rev(order.row) %>% str_replace('disease', 'dis.'))
df.plot <- df.plot %>% mutate(padj = replace(padj, padj > lim.lpa.lower, 1))
ggplot(df.plot, aes(x= cell, y=var, size=-log10(padj), color=Estimate)) + 
  geom_point() + 
  # scale_color_gradient2(low = "blue", mid = "white",  high = "red", space = "Lab", limit = c(-2,2)) +
  # scale_fill_brewer(palette = "RdBu", limit = c(-2,2)) +
  # scale_color_gradient(palette = "RdBu", space = "Lab", limit = c(-2,2)) +
  scale_color_gradientn(colours = rev(brewer.pal(10, "RdBu")), limit = c(lim.nest,lim.est)) +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  scale_size(limits = c(-log10(lim.lpa.lower),100)) +
  theme_bw()
ggsave(paste0(prefix.output, '_queryL2_glm_withproj.pdf'), width = 6.6, height = 5.2)


write_csv(df.res, paste0(prefix.output, '_queryL2_glm_withproj.csv'))
write_csv(df.output.wide, paste0(prefix.output, '_queryL2_output.csv'))

## cluster L1
list.df <- list()

for (project.name in project.name.all) {
  df.metadata <- read_csv(paste0('data/', project.name,'/', project.name,'_metadata.csv'))
  df.output <- read_csv(paste0('output/', project.name, '/',  project.name, '_queryL1_Reference_Mapping.csv'))
  
  df.output <- df.output %>% 
    left_join(df.metadata, by = "sample") %>%
    group_by(sample) %>%
    mutate(freq = n / sum(n), project = project.name,
           sample = paste(project.name, sample, sep = ":"))
  
  list.df[[project.name]] <- df.output
}

df.output <- rbindlist(list.df, fill = TRUE) %>% 
  as_tibble() %>%
  mutate(sex = str_to_lower(sex), age = age / 25)

df.output <- df.output %>%
  filter(!disease %in% blacklist_disease) %>%
  filter(!sample %in% blacklist_sample)

df.output$disease <- str_replace(df.output$disease, "COVID-19", "Covid")
# 230112 CD&SpA to CD and droped SpA
df.output$disease <- str_replace(df.output$disease, "CD&SpA", "CD")
df.output <- df.output %>% filter(disease != "SpA") %>%
  mutate(disease = relevel(factor(disease), ref = "HC"))
write_csv(df.output, paste0(prefix.output, '_tidy_queryL1_output.csv'))

n_total <- df.output %>% group_by(project) %>% summarise(n_total = sum(n))
n_hc <- df.output %>% group_by(project) %>% filter(disease == "HC") %>% summarise(n_HC = n_distinct(sample))
n_dis <- df.output %>% group_by(project) %>% filter(disease != "HC") %>% summarise(n_dis = n_distinct(sample))
df.stats <- n_total %>% left_join(n_hc, by = "project") %>% left_join(n_dis, by = "project")
write_csv(df.stats, paste0(prefix.output, '_clusterL1_stats.csv'))

p <- ggplot(df.output, aes(x=clusterL1, y=freq, color=disease), scale = "width") + 
  geom_quasirandom(dodge.width=1, size=0.5) +
  theme(axis.text.x = element_text(angle = 90))
p
ggsave(paste0(prefix.output, '_queryL1_ggbeeswarm.pdf'), width = 10, height = 4)

df.ntotal <- df.output %>% group_by(sample) %>% summarise(n_total = sum(n))
cols <- c("sample", "disease", "disease_duration", "condition", "age", "race", 
          "sex", "immunosuppressant", "immunosuppressant_duration", "project")
df.metadata <- df.output[,cols] %>% distinct()


df.output.wide <- df.output %>% 
  pivot_wider(values_from = n, names_from = clusterL1,
              id_cols = c(sample, disease, disease_duration, condition, age, race, sex,
                          immunosuppressant, immunosuppressant_duration, project), values_fill = 0) %>%
  left_join(df.ntotal, by='sample') %>%
  drop_na(age)

df.output.wide <- as.data.frame(df.output.wide)
df.output.wide <- within(df.output.wide, disease <- relevel(factor(disease), ref = "HC"))
df.output.wide <- within(df.output.wide, sex <- relevel(factor(sex), ref = "male"))

cells <- c("Tnaive", "Tcm", 'Tem', "Temra", "Treg")

list.res <- list()
for (cell in cells) {
  res <- glm(as.formula(paste0("cbind(", cell, ", n_total-", cell,") ~ disease + age + sex + project")),
             data = df.output.wide, family = binomial)
  # res <- glm(as.formula(paste0("cbind(", cell, ", n_total-", cell,") ~ disease + age + sex")),
  #            data = df.output.wide, family = binomial)
  list.res[[cell]] <- as_tibble(summary(res)$coefficients, rownames = "var") %>%
    mutate(cell = cell)
}

df.res <- rbindlist(list.res, fill = TRUE) %>% 
  as_tibble() %>%
  mutate(padj = p.adjust(`Pr(>|z|)`))

df.plot <- df.res %>%
  filter(str_detect(var, "^disease|^age|^sex")) %>%
  mutate(cell = factor(cell, levels = cells))

lim.est <- 1.8
lim.nest <- -1.8
lim.lpa <- 10^-100
lim.lpa.lower <- 0.05

df.plot <- df.plot %>% 
  mutate(Estimate=if_else(Estimate>lim.est, lim.est , Estimate)) %>% 
  mutate(Estimate=if_else(Estimate< lim.nest, -lim.est, Estimate)) %>%
  mutate(padj=if_else(padj<lim.lpa, lim.lpa, padj)) %>%
  mutate(var=fct_relevel(var, "sexfemale", "age"))

order.row <- df.plot %>% group_by(var) %>% filter(`z value` == max(`z value`)) %>% arrange(cell) %>% select(var)
order.row <- as.character(order.row$var)
order.row <- order.row[!order.row %in% c("age", "sexfemale")]
order.row <- append(order.row, c('age', 'sexfemale'))
# df.plot$var <- factor(df.plot$var, levels = rev(order.row))
df.plot$var <- str_replace(df.plot$var, 'disease', 'dis.')
df.plot$var <- factor(df.plot$var, levels = rev(order.row) %>% str_replace('disease', 'dis.'))
df.plot <- df.plot %>% mutate(padj = replace(padj, padj > lim.lpa.lower, 1))
ggplot(df.plot, aes(x= cell, y=var, size=-log10(padj), color=Estimate)) + 
  geom_point() + 
  # scale_color_gradient2(low = "blue", mid = "white",  high = "red", space = "Lab", limit = c(-2,2)) +
  scale_color_gradientn(colours = rev(brewer.pal(10, "RdBu")), limit = c(lim.nest,lim.est)) +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  scale_size(limits = c(-log10(lim.lpa.lower),100)) +
  theme_bw()
ggsave(paste0(prefix.output, '_queryL1_glm_withproj.pdf'), width = 4.5, height = 5.2)


write_csv(df.res, paste0(prefix.output, '_queryL1_glm_withproj.csv'))
write_csv(df.output.wide, paste0(prefix.output, '_queryL1_output.csv'))


### NMF
cells <- c("Tnaive", "TnaiveAct", "TnaiveMX1", "TnaiveSOX4","TcmTh0","TcmTh0Act", 
           "TcmTfh", "TcmPHLDA3", "TcmTh17", "TcmTh2",
           "TemTh1pre", "TemTh1", "TemTh117", "TemTph", "TemraTh1",
           "TregNaive", "TregAct", "TregEff")

list.df <- list()

for (project.name in project.name.all) {
  df.metadata <- read_csv(paste0('data/', project.name,'/', project.name,'_metadata.csv'))
  df.output <- read_csv(paste0('output/', project.name, '/',  project.name, '_CD4T_NMF_arranged.csv'))
  
  df.output <- df.output %>% 
    left_join(df.metadata, by = "sample", suffix = c(".x", "")) %>%
    mutate(project = project.name,
           sample = paste(project.name, sample, sep = ":")) %>%
    filter(!disease %in% blacklist_disease) %>%
    filter(!sample %in% blacklist_sample)
  
  list.df[[project.name]] <- df.output
}

df.output <- rbindlist(list.df, fill = TRUE) %>% 
  as_tibble() %>% 
  drop_na(disease) %>% 
  mutate(sex = str_to_lower(sex), age = age / 25) %>%
  mutate(disease = relevel(factor(disease), ref = "HC")) %>%
  mutate(sex = relevel(factor(sex), ref = "male"))
df.output$disease <- str_replace(df.output$disease, "COVID-19", "Covid")
df.output$disease <- str_replace(df.output$disease, "CD&SpA", "CD")
df.output <- df.output %>% filter(disease != "SpA") %>%
  mutate(disease = relevel(factor(disease), ref = "HC"))

df.output <- within(df.output, disease <- relevel(factor(disease), ref = "HC"))
# df.output <- within(df.output, sex <- relevel(factor(sex), ref = "male"))
write_csv(df.output, paste0(prefix.output, '_tidy_NMFL2_output.csv'))

df.output %>% filter(clusterL2=="TnaiveMX1") %>% ggplot(aes(NMF_1)) + geom_histogram()

# cases
df.output %>% filter(disease != "HC") %>% select(sample) %>% distinct() %>% count()

# controls
df.output %>% filter(disease == "HC") %>% select(sample) %>% distinct() %>% count()

# num of cells
df.output %>% count()


# save matrix
write_csv(df.output, paste0(prefix.output, '_queryL2_output_withNMF.csv'))

df.plot <- df.output %>% group_by(clusterL1, disease) %>%
  summarise(NMF_0=mean(NMF_0), NMF_1=mean(NMF_1), NMF_2=mean(NMF_2),
            NMF_3=mean(NMF_3), NMF_4=mean(NMF_4), NMF_5=mean(NMF_5),
            NMF_6=mean(NMF_6), NMF_7=mean(NMF_7), NMF_8=mean(NMF_8),
            NMF_9=mean(NMF_9), NMF_10=mean(NMF_10), NMF_11=mean(NMF_11)) %>%
  mutate(description = paste0(clusterL1,'_',disease))

col.labels <- c('NMF0 Cytotoxic-F', 'NMF1 Treg-F', 'NMF2 Th17-F', 'NMF3 Naive-F', 
                'NMF4 Act-F', 'NMF5 TregEff/Th2-F', 'NMF6 Tfh-F', 'NMF7 IFN-F', 'NMF8 Cent.Mem.-F',
                'NMF9 Thy.Emi.-F', 'NMF10 Tissue-F', 'NMF11 Th1-F')
col.factors <- c('NMF_0', 'NMF_1', 'NMF_2',
                 'NMF_3', 'NMF_4', 'NMF_5',
                 'NMF_6', 'NMF_7', 'NMF_8',
                 'NMF_9', 'NMF_10', 'NMF_11')

p <- pheatmap(df.plot[col.factors],
              # color = cividis(100),
              cluster_cols = FALSE, cluster_rows = FALSE,
              labels_row = df.plot$description, labels_col = col.labels,
              scale = "column")
p

pdf(paste0(prefix.output, '_NMF_clusterL1_heatmap.pdf'), width = 6, height = 10)
p
dev.off()

df.plot <- df.output %>% group_by(clusterL2, disease) %>%
  summarise(NMF_0=mean(NMF_0), NMF_1=mean(NMF_1), NMF_2=mean(NMF_2),
            NMF_3=mean(NMF_3), NMF_4=mean(NMF_4), NMF_5=mean(NMF_5),
            NMF_6=mean(NMF_6), NMF_7=mean(NMF_7), NMF_8=mean(NMF_8),
            NMF_9=mean(NMF_9), NMF_10=mean(NMF_10), NMF_11=mean(NMF_11)) %>%
  mutate(description = paste0(clusterL2,'_',disease))

p <- pheatmap(df.plot[c('NMF_0', 'NMF_1', 'NMF_2',
                        'NMF_3', 'NMF_4', 'NMF_5',
                        'NMF_6', 'NMF_7', 'NMF_8',
                        'NMF_9', 'NMF_10', 'NMF_11')],
              # color = cividis(100),
              cluster_cols = FALSE, cluster_rows = FALSE,
              labels_row = df.plot$description, labels_col = col.labels,
              scale = "column")
p

pdf(paste0(prefix.output, '_NMF_clusterL2_heatmap.pdf'), width = 6, height = 32)
p
dev.off()
# 
# list.res <- list()
# for (cell in cells) {
#   for (comp in colnames(df.output) %>% str_subset("NMF_")){
#     # d <- df.output %>% filter(clusterL2 == cell) %>% 
#     #   group_by(sample) %>% 
#     #   summarise(NMF = mean(.data[[comp]]), disease=first(disease), 
#     #             age=first(age), sex=first(sex), project=first(project))
#     # res <- glm(as.formula("NMF ~ disease + age + sex + project"),
#     #            data = d, family = Gamma)
#     res <- glm(as.formula(paste0(comp, " ~ disease + age + sex")),
#                data = df.output %>% filter(clusterL2 == cell))
#     list.res[[paste(cell, comp)]] <- as_tibble(summary(res)$coefficients, rownames = "var") %>%
#       mutate(cell = cell, component = comp)
#   }
# }
# 
# df.res <- rbindlist(list.res, fill = TRUE) %>% 
#   as_tibble() %>%
#   mutate(padj = p.adjust(`Pr(>|t|)`))
# 
# 
# for (c in cells) {
#   df.plot <- df.res %>%
#     filter(cell == c) %>%
#     filter(str_detect(var, "^disease|^age|^sex")) %>%
#     filter(padj < 0.1) %>% 
#     mutate(cell = factor(cell, levels = cells))
#   
#   lim.est <- 0.5
#   lim.nest <- -0.5
#   lim.lpa <- 10^-100
#   
#   df.plot <- df.plot %>% 
#     mutate(Estimate=if_else(Estimate>lim.est, lim.est , Estimate)) %>% 
#     mutate(Estimate=if_else(Estimate< lim.nest, -lim.est, Estimate)) %>%
#     mutate(padj=if_else(padj<lim.lpa, lim.lpa, padj))
#   
#   if (dim(df.plot)[1] > 0){
#     ggplot(df.plot, aes(x= component, y=var, size=-log10(padj), color=Estimate)) + 
#       geom_point() + 
#       scale_color_gradient2(low = "blue", mid = "white",  high = "red", space = "Lab", limit = c(-0.5,0.5)) +
#       scale_x_discrete(guide = guide_axis(angle = 90)) +
#       theme_classic()
#     ggsave(paste0(prefix.output, '_NMF_queryL2', c, '_glm_withoutproj.pdf'), width = 6, height = 5)
#   }
# }
# 
# write_csv(df.res, paste0(prefix.output, '_NMF_queryL2_glm_withoutproj.csv'))


list.res <- list()
for (cell in cells) {
  for (comp in colnames(df.output) %>% str_subset("NMF_")){
    # d <- df.output %>% filter(clusterL2 == cell) %>% 
    #   group_by(sample) %>% 
    #   summarise(NMF = mean(.data[[comp]]), disease=first(disease), 
    #             age=first(age), sex=first(sex), project=first(project))
    # res <- glm(as.formula("NMF ~ disease + age + sex + project"),
    #            data = d, family = Gamma)
    res <- glm(as.formula(paste0(comp, " ~ disease + age + sex + project")),
               data = df.output %>% filter(clusterL2 == cell))
    list.res[[paste(cell, comp)]] <- as_tibble(summary(res)$coefficients, rownames = "var") %>%
      mutate(cell = cell, component = comp)
  }
}

df.res <- rbindlist(list.res, fill = TRUE) %>% 
  as_tibble() %>%
  mutate(padj = p.adjust(`Pr(>|t|)`))


lim.est <- 0.3
lim.nest <- -0.3
lim.lpa <- 10^-100
lim.scaledEstimate <- 3

for (c in cells) {
  df.plot <- df.res %>%
    filter(cell == c) %>%
    filter(str_detect(var, "^disease|^age|^sex")) %>%
    mutate(cell = factor(cell, levels = cells))
  
  df.plot <- df.plot %>% 
    mutate(Estimate=if_else(Estimate>lim.est, lim.est , Estimate)) %>% 
    mutate(Estimate=if_else(Estimate< lim.nest, -lim.est, Estimate)) %>%
    mutate(padj=if_else(padj<lim.lpa, lim.lpa, padj)) %>% 
    mutate(padj=if_else(padj>0.1, 0, padj))
  
  df.plot <- df.plot %>% mutate(component=as_factor(component)) %>% 
    mutate(component=fct_relevel(component, col.factors)) %>%
    mutate(component=fct_recode(component, 
                                'NMF0 Cytotoxic-F'='NMF_0', 'NMF1 Treg-F'='NMF_1', 
                                'NMF2 Th17-F'='NMF_2', 'NMF3 Naive-F'='NMF_3', 
                                'NMF4 Act-F'='NMF_4', 'NMF5 TregEff/Th2-F'='NMF_5', 
                                'NMF6 Tfh-F'='NMF_6', 'NMF7 IFN-F'='NMF_7', 
                                'NMF8 Cent.Mem.-F'='NMF_8', 
                                'NMF9 Thy.Emi.-F'='NMF_9', 
                                'NMF10 Tissue-F'='NMF_10',
                                'NMF11 Th1-F'='NMF_11')) %>%
    mutate(var=fct_relevel(var, "sexfemale", "age"))
  

  df.plot$var <- str_replace(df.plot$var, 'disease', 'dis.')
  df.plot$var <- factor(df.plot$var, levels = rev(order.row) %>% str_replace('disease', 'dis.'))
  
  df.plot.int <- df.res %>%
    filter(var == "(Intercept)") %>%
    mutate(component=as_factor(component)) %>% 
    mutate(component=fct_relevel(component, col.factors)) %>%
    mutate(component=fct_recode(component, 
                                'NMF0 Cytotoxic-F'='NMF_0', 'NMF1 Treg-F'='NMF_1', 
                                'NMF2 Th17-F'='NMF_2', 'NMF3 Naive-F'='NMF_3', 
                                'NMF4 Act-F'='NMF_4', 'NMF5 TregEff/Th2-F'='NMF_5', 
                                'NMF6 Tfh-F'='NMF_6', 'NMF7 IFN-F'='NMF_7', 
                                'NMF8 Cent.Mem.-F'='NMF_8', 
                                'NMF9 Thy.Emi.-F'='NMF_9', 
                                'NMF10 Tissue-F'='NMF_10',
                                'NMF11 Th1-F'='NMF_11')) %>%
    group_by(component) %>%
    mutate(z_scaled_Estimate = as.numeric(scale(Estimate))) %>%
    ungroup() %>%
    mutate(z_scaled_Estimate=if_else(z_scaled_Estimate>lim.scaledEstimate, lim.scaledEstimate, z_scaled_Estimate)) %>% 
    # mutate(z_scaled_Estimate=if_else(z_scaled_Estimate<-lim.scaledEstimate, -lim.scaledEstimate, z_scaled_Estimate)) %>% 
    filter(cell == c) %>%
    mutate(cell = factor(cell, levels = cells))
  
  if (dim(df.plot)[1] > 0){
    # Create the dotplot
    dotplot <- ggplot(df.plot, aes(x= component, y=var, size=-log10(padj), color=Estimate)) + 
      geom_point() +
      scale_color_gradientn(colours = rev(brewer.pal(10, "RdBu")), limit = c(lim.nest,lim.est)) +
      scale_x_discrete(guide = guide_axis(angle = 90)) +
      scale_size(range = c(0, 4), limits = c(0,100)) + 
      xlab("") + ylab("") +
      labs(color = "Coef.") +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5, colour = "black"),
            axis.text = element_text(color = "black"),
            axis.title = element_text(color = "black"),
            legend.text = element_text(color = "black"),
            legend.title = element_text(color = "black"))
    
    # # Save the dotplot
    # ggsave(paste0(prefix.output, '_NMF_queryL2', c, '_glm_withproj.pdf'), width = 6, height = 5)
    
    # Create the heatmap
    heatmap <- ggplot(df.plot.int, aes(x = component, y = var, fill = z_scaled_Estimate)) +
      geom_tile() +
      ggtitle(c) +
      scale_fill_gradientn(colours = rev(brewer.pal(10, "RdGy")), limit = c(-lim.scaledEstimate,lim.scaledEstimate)) +
      theme_minimal() +
      labs(fill = "Intercept") +
      theme(
        axis.title.x = element_blank(), # Remove x-axis title
        axis.title.y = element_blank(), # Remove y-axis title
        axis.text.x = element_blank(),  # Remove x-axis text
        axis.text.y = element_blank(),  # Remove y-axis text
        panel.grid.major = element_blank(), # Remove major grid
        panel.grid.minor = element_blank(), # Remove minor grid
        legend.position = "right",
        legend.text = element_text(color = "black"),
        legend.title = element_text(color = "black")
      )
    
    # Combine the plots
    heatmap / dotplot + plot_layout(heights = c(1, 8))
    ggsave(paste0(prefix.output, '_NMF_queryL2', c, '_glm_withproj_withscaledintercept.pdf'), width = 6, height = 6)
  }
}

 write_csv(df.res, paste0(prefix.output, '_NMF_queryL2_glm_withproj.csv'))
　

