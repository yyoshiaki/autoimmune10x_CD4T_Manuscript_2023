## ----library and source---------------------------------------------------------------------------
library(symphony)
library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(tidyverse)
library(Matrix)
library(sctransform)

source("./script/symphony/vignettes/utils_seurat.R")
source('./script/functions.R')


## ----parameter setting----------------------------------------------------------------------------
project.name <- c("GSE194315")
prefix_data <- c("./data/")
prefix <- paste0("./output/", project.name, "/", project.name)


## ----load h5seurat--------------------------------------------------------------------------------
load(paste0(prefix_data, "ref_Reference_Mapping.RData"))


## ----Load data for query--------------------------------------------------------------------------
query_obj <- readRDS(paste0(prefix, "_CD4T_AssayData.rds"))


## ----process--------------------------------------------------------------------------------------
reference_mapping(ref, query_obj, prefix)


## ----save.image-----------------------------------------------------------------------------------
prefix = paste0("./output/", project.name, "/")
save.image(file = paste0(prefix, "220423_reference_mapping.RData"))


## ----sessionInfo----------------------------------------------------------------------------------
sessionInfo()
