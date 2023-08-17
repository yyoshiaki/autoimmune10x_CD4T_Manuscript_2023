## ----Ensure, eval=FALSE, include=FALSE------------------------------------------------------------
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


## ----library and source---------------------------------------------------------------------------
library(Seurat)
library(SeuratData)
library(Azimuth)
library(patchwork)
library(tidyverse)
library(sctransform)

setwd('/home/rstudio/autoimmune_10x')
source('./script/functions.R')


## ----parameter setting----------------------------------------------------------------------------
project.name <- c("3993994")
prefix_data <- paste0("./data/", project.name, "/")
prefix <- paste0("./output/", project.name, "/", project.name)


## ----Load_the_reference---------------------------------------------------------------------------
reference <- LoadReference(path = "https://seurat.nygenome.org/azimuth/references/v1.0.0/human_pbmc")


## ----reading data---------------------------------------------------------------------------------
query.data <- readRDS(paste0(prefix_data, "counts.rds"))
query <- CreateSeuratObject(counts = query.data,
                            project = project.name,
                            min.cells = 3,
                            min.features = 200)


## ----loading batch_file and metadata--------------------------------------------------------------
metadata <- read_tsv(paste0(prefix_data, "metadata.txt"))
#meta <- metadata %>%
#  select(sample, sampleType) %>%
#  group_by(sample) %>%
#  distinct()
#write.csv(meta,
#          file = paste0(prefix_data, "metadata.csv"),
#          row.names = FALSE,
#          quote = FALSE)

metadata <- metadata %>%
  select(cell, sample)
m <- tibble(query@assays[["RNA"]]@counts@Dimnames[[2]])
colnames(m) <- c("cell")
metadata <- metadata %>%
  right_join(m, by = c("cell"))
rownames(metadata) <- metadata$cell

query <- AddMetaData(query,
                     metadata = metadata)

meta <- read_csv(paste0(prefix_data, project.name, "_metadata.csv"))
metadata <- metadata %>%
  group_by(sample) %>%
  right_join(meta, by = c("sample"))
metadata <- metadata %>%
  relocate(cell,
           sample,
           disease,
           disease_duration,
           condition,
           age,
           race,
           sex,
           immunosuppressant,
           immunosuppressant_duration,
           batch)
rownames(metadata) <- metadata$cell

query <- query %>%
  subset(sample %in% meta$sample)
query <- AddMetaData(query,
                     metadata = metadata)


## ----process--------------------------------------------------------------------------------------
extract_cells(query, reference, prefix)


## ----save.image-----------------------------------------------------------------------------------
prefix = paste0("./output/", project.name, "/")
save.image(file = paste0(prefix, "220423_cell_extraction_Azimuth.RData"))


## ----sessionInfo----------------------------------------------------------------------------------
sessionInfo()

