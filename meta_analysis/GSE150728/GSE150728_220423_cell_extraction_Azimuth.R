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
library(Matrix)

setwd('/home/rstudio/autoimmune_10x')
source('./script/functions.R')

## ----parameter setting----------------------------------------------------------------------------------------
project.name <- "GSE150728"
prefix_data <- paste0("./data/", project.name, "/")
#prefix <- paste0("./output/", project.name, "/", project.name)


## ----Load_the_reference---------------------------------------------------------------------------------------
reference <- LoadReference(path = "https://seurat.nygenome.org/azimuth/references/v1.0.0/human_pbmc")


## ----reading data---------------------------------------------------------------------------------------------
## ----process--------------------------------------------------------------------------------------------------
filename <- read_csv(paste0(prefix_data, "GSE150728_filename.csv"))
metadata <- read_csv(paste0(prefix_data, "GSE150728_metadata.csv"))

for (fn in filename$filename) {
  print(fn)
  f <- filename %>%
    subset(filename == fn)
  m <- metadata %>%
    subset(sample == f$sample)
  
  so <- readRDS(paste0(prefix_data, fn))
  so <- CreateSeuratObject(counts = so[["exon"]],
                           project = project.name,
                           min.cells = 3,
                           min.features = 200)
  cells <- so@assays[["RNA"]]@counts@Dimnames[[2]] %>%
    tibble()
  colnames(cells) <- c("cells")
  cells <- cells %>%
    mutate(sample = f$sample)
  cells <- cells %>%
    group_by(sample) %>%
    left_join(m, by = "sample")
  rownames(cells) <- cells$cells
  so <- AddMetaData(so,
                    metadata = cells)
  
  prefix = paste0("./output/", project.name, "/", f$sample)
  
  extract_cells(so, reference, prefix)
  }


## ----save.image-----------------------------------------------------------------------------------------------
prefix = paste0("./output/", project.name, "/")
save.image(file = paste0(prefix, "220423_CD4Tcell_extraction_Azimuth.RData"))


## ----sessionInfo----------------------------------------------------------------------------------------------
sessionInfo()

