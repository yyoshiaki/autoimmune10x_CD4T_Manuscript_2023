## ----Ensure, eval=FALSE, include=FALSE------------------------------------------------------------------------
## # Ensure Seurat v4.0 or higher is installed
## if (packageVersion(pkg = "Seurat") < package_version(x = "4.0.0")) {
##   stop("Mapping datasets requires Seurat v4 or higher.", call. = FALSE)
## }
## 
## # Ensure glmGamPoi is installed
## if (!requireNamespace("glmGamPoi", quietly = TRUE)) {
##   if (!requireNamespace("BiocManager", quietly = TRUE)) {
##     BiocManager::install("glmGamPoi")
##   }
## }
## 
## # Ensure Azimuth is installed
## if (packageVersion(pkg = "Azimuth") < package_version(x = "0.3.1")) {
##   stop("Please install azimuth - remotes::install_github('satijalab/azimuth')", call. = FALSE)
## }


## ----library--------------------------------------------------------------------------------------------------
library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(Azimuth)
library(patchwork)
library(tidyverse)
library(sctransform)

setwd('/home/rstudio/autoimmune_10x')
source('./script/functions.R')

## ----parameter setting----------------------------------------------------------------------------------------
project.name <- "577860"
prefix_data <- paste0("./data/", project.name, "/")
#prefix <- paste0("./output/", project.name, "/", project.name)

## ----Load_the_reference---------------------------------------------------------------------------------------
reference <- LoadReference(path = "https://seurat.nygenome.org/azimuth/references/v1.0.0/human_pbmc")


## ----reading data---------------------------------------------------------------------------------------------
## ----process--------------------------------------------------------------------------------------------------
sample <- read_csv(paste0(prefix_data, "577860_metadata.csv"))

for (s in sample$sample) {
  print(s)
  prefix_data <- paste0("./data/", project.name, "/raw_umi/", s, "_PBMC")
  
  query.data <- Read10X(data.dir = prefix_data)
  
  prefix = paste0("./output/", project.name, "/", s)
  
  q <- CreateSeuratObject(counts = query.data,
                          project = project.name,
                          assay = "RNA",
                          min.cells = 3,
                          min.features = 200)
  
  meta <- filter(sample, sample == s)
  rowname <- tibble(q@assays[["RNA"]]@counts@Dimnames[[2]])
  rowname <- mutate(rowname, batch = s)
  rowname <- rowname %>%
    group_by(batch) %>%
    left_join(meta, by = c("batch"))

  rownames(rowname) <- rowname$`q@assays[["RNA"]]@counts@Dimnames[[2]]`
  
  q <- AddMetaData(q, rowname)
  
  extract_cells(q, reference, prefix)
  }


## ----save.image-----------------------------------------------------------------------------------------------
prefix = paste0("./output/", project.name, "/")
save.image(file = paste0(prefix, "220423_CD4Tcell_extraction_Azimuth.RData"))


## ----sessionInfo----------------------------------------------------------------------------------------------
sessionInfo()

