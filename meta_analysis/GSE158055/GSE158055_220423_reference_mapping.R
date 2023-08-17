## ----library and source---------------------------------------------------------------------------
library(symphony)
library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(tidyverse)
library(Matrix)
library(sctransform)

setwd('/home/rstudio/autoimmune_10x')
source("./script/symphony/vignettes/utils_seurat.R")
source('./script/functions.R')


## ----parameter setting----------------------------------------------------------------------------
project.name <- c("GSE158055")
prefix_data <- c("./data/")
prefix <- paste0("./output/", project.name, "/", project.name)


## ----load h5seurat--------------------------------------------------------------------------------
load(paste0(prefix_data, "ref_Reference_Mapping.RData"))


## ----Load data for query--------------------------------------------------------------------------
## ----process--------------------------------------------------------------------------------------
sample <- read_csv(paste0(prefix, "_processed_sample_cell_extraction.csv"))

ss <- 0
for (i in sample$sample) {
  print(i)
  prefix = paste0("./output/", project.name, "/", i)
  if (file.exists(paste0(prefix, "_CD4T_AssayData.rds"))) {
    query_obj <- readRDS(paste0(prefix, "_CD4T_AssayData.rds"))
    if (query_obj@Dim[2] >= 500) {
      reference_mapping(ref, query_obj, prefix)
      if (ss==0) {
        ii <- tibble(i)
        colnames(ii) <- "sample"
        ss <- ii
      } else {
        ii <- tibble(i)
        colnames(ii) <- "sample"
        ss <- ss %>% bind_rows(ii)
      }
    } else {
      print(paste0(i, " was removed from processing because cell number is <500"))
    }
  } else {
    print("skipped...")
  }
}
write_csv(ss, file = paste0("./output/", project.name, "/", project.name, "_processed_sample_reference_mapping.csv"))

df <- 0
for (s in ss$sample) {
  prefix = paste0("./output/", project.name, "/", s)
  f <- paste0(prefix, "_queryL1_Reference_Mapping.csv")
  if (file.exists(f)){
    # print(f)
    d <- read.csv(f)
    if (df==0) {
      df <- d
    } else {
      df <- df %>% bind_rows(d)
    }
  }
}
write_csv(df, file = paste0("./output/", project.name, "/", project.name, "_queryL1_Reference_Mapping.csv"))

df <- 0
for (s in ss$sample) {
  prefix = paste0("./output/", project.name, "/", s)
  f <- paste0(prefix, "_queryL2_Reference_Mapping.csv")
  if (file.exists(f)){
    # print(f)
    d <- read.csv(f)
    if (df==0) {
      df <- d
    } else {
      df <- df %>% bind_rows(d)
    }
  }
}
write_csv(df, file = paste0("./output/", project.name, "/", project.name, "_queryL2_Reference_Mapping.csv"))

meta <- read_csv(paste0(prefix_data, project.name, "/",project.name, "_metadata.csv"))
meta_processed <- meta %>%
  right_join(ss, by = c("sample"))
write_csv(meta_processed, file = paste0("./output/", project.name, "/", project.name, "_metadata_processed_sample_reference_mapping.csv"))


## ----save.image-----------------------------------------------------------------------------------
#prefix = paste0("./output/", project.name, "/")
#save.image(file = paste0(prefix, "220423_reference_mapping.RData"))


## ----sessionInfo----------------------------------------------------------------------------------
sessionInfo()

