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
project.name <- c("GSE157278")

prefix_data <- paste0("./data/", project.name, "/")

prefix <- paste0("./output/", project.name, "/", project.name)


## ----Load_the_reference---------------------------------------------------------------------------
reference <- LoadReference(path = "https://seurat.nygenome.org/azimuth/references/v1.0.0/human_pbmc")


## ----reading data---------------------------------------------------------------------------------
query.data <- Read10X(data.dir = prefix_data)
query <- CreateSeuratObject(counts = query.data,
                            project = project.name,
                            min.cells = 3,
                            min.features = 200)
query


## ----pbmc@assays----------------------------------------------------------------------------------
query@assays[["RNA"]]@counts[1:10,1:20]


## ----loading batch_file and metadata--------------------------------------------------------------
batch <- read_table(paste0(prefix_data, "/cell_batch.tsv.gz"))

batch_select <- query@assays[["RNA"]]@data@Dimnames[[2]] %>%
  tibble()
colnames(batch_select) <- c("Cell")
batch <- right_join(batch,
                    batch_select,
                    by = "Cell")

metadata <- read.csv(paste0(prefix_data, project.name, "_metadata.csv"))

batch <- batch %>%
  group_by(batch) %>%
  left_join(metadata,
            by = c("batch"))

batch <- batch %>%
  relocate(Cell,
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

rownames(batch) <- batch$Cell

query <- AddMetaData(query,
                     metadata = batch)


## ----process--------------------------------------------------------------------------------------
extract_cells(query, reference, prefix)


## ----save.image-----------------------------------------------------------------------------------
prefix = paste0("./output/", project.name, "/")
save.image(file = paste0(prefix, "220423_cell_extraction_Azimuth.RData"))


## ----sessionInfo----------------------------------------------------------------------------------
sessionInfo()

