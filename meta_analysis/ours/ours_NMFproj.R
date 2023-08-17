## ----library and source---------------------------------------------------------------------------
library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(tidyverse)
library(Matrix)
library(sctransform)
library(pheatmap)
library(viridis)

setwd('~/autoimmune_10x/')
# source("./script/symphony/vignettes/utils_seurat.R")
# source('./script/functions.R')


## ----parameter setting----------------------------------------------------------------------------
project.name <- c("ours")
prefix_data <- c("./data/")
prefix <- paste0("./output/", project.name, "/", project.name)

## ----NMF-----------------------------------------------------------------------------------------

load(paste0(prefix_data, "ref_Reference_Mapping.RData"))
# ref <- LoadH5Seurat("./data/res.forseurat.v4.h5seurat", assays = "data")

cmd <- paste("/home/rstudio/anaconda3/bin/python script/NMFprojection/NMFprojection.py", 
             "./data/ours/ours_CD4T_X.tsv",
             "--outputprefix", paste0(prefix, '_CD4T'),
             "script/NMFprojection/data/NMF.W.CD4T.csv.gz")

# system(cmd)

# L1 
# load(paste0("./output/", project.name, "/", project.name, "_queryL1_Reference_Mapping.RData"))
# load(paste0("./output/", project.name, "/", project.name, "_queryL2_Reference_Mapping.RData"))

df.proj <- read.csv(paste0(prefix, '_CD4T_projection.csv'), row.names = 1, check.names = FALSE)
df.proj <- t(df.proj)
head(df.proj)

df <- ref$meta_data %>% 
  mutate(condition = disease) %>%
  select(c(sample, disease, condition, age, sex, clusterL1,  clusterL2)) 
df <- merge(df, df.proj, by = 'row.names', all = TRUE)
df <- df[-c(1)]
df <- df %>% drop_na(sample)
head(df)

write_csv(df, paste0(prefix, '_CD4T_NMF_arranged.csv'))

df.plot <- df %>% group_by(clusterL1, disease) %>%
  summarise(NMF_0=mean(NMF_0), NMF_1=mean(NMF_1), NMF_2=mean(NMF_2),
            NMF_3=mean(NMF_3), NMF_4=mean(NMF_4), NMF_5=mean(NMF_5),
            NMF_6=mean(NMF_6), NMF_7=mean(NMF_7), NMF_8=mean(NMF_8),
            NMF_9=mean(NMF_9), NMF_10=mean(NMF_10), NMF_11=mean(NMF_11)) %>%
  mutate(description = paste0(clusterL1,'_',disease))

col.labels <- c('NMF0 Cytotoxic', 'NMF1 Treg', 'NMF2 Th17', 'NMF3 Naiveness', 
                'NMF4 Act', 'NMF5 Th2', 'NMF6 Tfh', 'NMF7 IFN', 'NMF8 Cent. Mem.',
                'NMF9 Thymic Emi.', 'NMF10 Tissue', 'NMF11 Th1')
p <- pheatmap(df.plot[c('NMF_0', 'NMF_1', 'NMF_2',
                        'NMF_3', 'NMF_4', 'NMF_5',
                        'NMF_6', 'NMF_7', 'NMF_8',
                        'NMF_9', 'NMF_10', 'NMF_11')],
              # color = cividis(100),
              cluster_cols = FALSE, cluster_rows = FALSE,
              labels_row = df.plot$description, labels_col = col.labels,
              scale = "column")
p

pdf(paste0(prefix, '_CD4T_NMF_clusterL1_heatmap.pdf'), width = 6, height = 4)
p
dev.off()

df.plot <- df %>% group_by(clusterL2, disease) %>%
  summarise(NMF_0=mean(NMF_0), NMF_1=mean(NMF_1), NMF_2=mean(NMF_2),
            NMF_3=mean(NMF_3), NMF_4=mean(NMF_4), NMF_5=mean(NMF_5),
            NMF_6=mean(NMF_6), NMF_7=mean(NMF_7), NMF_8=mean(NMF_8),
            NMF_9=mean(NMF_9), NMF_10=mean(NMF_10), NMF_11=mean(NMF_11)) %>%
  mutate(description = paste0(clusterL2,'_',disease))

col.labels <- c('NMF0 Cytotoxic', 'NMF1 Treg', 'NMF2 Th17', 'NMF3 Naiveness', 
                'NMF4 Act', 'NMF5 Th2', 'NMF6 Tfh', 'NMF7 IFN', 'NMF8 Cent. Mem.',
                'NMF9 Thymic Emi.', 'NMF10 Tissue', 'NMF11 Th1')
p <- pheatmap(df.plot[c('NMF_0', 'NMF_1', 'NMF_2',
                        'NMF_3', 'NMF_4', 'NMF_5',
                        'NMF_6', 'NMF_7', 'NMF_8',
                        'NMF_9', 'NMF_10', 'NMF_11')],
              # color = cividis(100),
              cluster_cols = FALSE, cluster_rows = FALSE,
              labels_row = df.plot$description, labels_col = col.labels,
              scale = "column")
p

pdf(paste0(prefix, '_CD4T_NMF_clusterL2_heatmap.pdf'), width = 6, height = 8)
p
dev.off()

