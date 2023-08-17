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
project.name <- "GSE150728"
prefix_data <- c("./data/")
#prefix <- paste0("./output/", project.name, "/", project.name)

df.metadata <- read_csv(paste0(prefix_data, project.name, '/', project.name, "_metadata.csv"))

## ----load h5seurat--------------------------------------------------------------------------------
load(paste0(prefix_data, "ref_Reference_Mapping.RData"))

## ----process--------------------------------------------------------------------------------------
for (sample in df.metadata$sample) {
  print(sample)
  prefix = paste0("./output/", project.name, "/", sample)
  if (file.exists(paste0(prefix, "_CD4T_AssayData.rds"))) {
    query_obj <- readRDS(paste0(prefix, "_CD4T_AssayData.rds"))
    reference_mapping(ref, query_obj, prefix)
  } else {
    print("skipped...")
  }
}

df <- 0
for (sample in df.metadata$sample) {
  prefix = paste0("./output/", project.name, "/", sample)
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
for (sample in df.metadata$sample) {
  prefix = paste0("./output/", project.name, "/", sample)
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

## ----save.image-----------------------------------------------------------------------------------
prefix = paste0("./output/", project.name, "/")
save.image(file = paste0(prefix, "220423_reference_mapping.RData"))


## ----sessionInfo----------------------------------------------------------------------------------
sessionInfo()

