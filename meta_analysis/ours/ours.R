## ----library and source---------------------------------------------------------------------------
library(Seurat)
library(tidyverse)

source("./script/symphony/vignettes/utils_seurat.R")
source('./script/functions.R')


## ----parameter setting----------------------------------------------------------------------------
project.name <- c("ours")
prefix_data <- c("./data/", project.name)
prefix <- paste0("./output/", project.name, "/", project.name)

load(paste0(prefix_data, "ref_Reference_Mapping.RData"))

df.metadata <- ref$meta_data

df.metadata <- df.metadata %>%
  mutate(disease_duration = NA,
         condition = NA,
         race = "EastAsian",
         immunosuppressant = 0,
         immunosuppressant_duration = 0,
         batch = ID) %>%
  select(sample,
           disease,
           disease_duration,
           condition,
           age,
           race,
           sex,
           immunosuppressant,
           immunosuppressant_duration,
           batch) %>%
  distinct()

write_csv(df.metadata, paste0(prefix_data, project.name, "_metadata.csv")) 


df.output <- ref$meta_data

queryL1.meta.data <- df.output %>% count(sample, clusterL1) %>% 
  mutate(clusterL1_prob.mean = 1,
         clusterL1_prob.max = 1,
         clusterL1_prob.min = 1) %>% 
  select(sample,clusterL1,n,clusterL1_prob.mean,clusterL1_prob.max,clusterL1_prob.min)

write.csv(queryL1.meta.data,
          file = paste0(prefix, "_queryL1_Reference_Mapping.csv"),
          row.names = FALSE,
          quote = FALSE)

queryL2.meta.data <- df.output %>% count(sample, clusterL2) %>% 
  mutate(clusterL2_prob.mean = 1,
         clusterL2_prob.max = 1,
         clusterL2_prob.min = 1) %>% 
  select(sample,clusterL2,n,clusterL2_prob.mean,clusterL2_prob.max,clusterL2_prob.min)

write.csv(queryL2.meta.data,
          file = paste0(prefix, "_queryL2_Reference_Mapping.csv"),
          row.names = FALSE,
          quote = FALSE)

## ----sessionInfo----------------------------------------------------------------------------------
sessionInfo()