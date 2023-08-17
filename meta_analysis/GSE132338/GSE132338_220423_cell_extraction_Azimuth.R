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
library(Matrix)

source('./script/functions.R')


## ----parameter setting----------------------------------------------------------------------------
project.name <- c("GSE132338")
prefix_data <- paste0("./data/", project.name, "/")
prefix <- paste0("./output/", project.name, "/", project.name)


## ----Load_the_reference---------------------------------------------------------------------------
reference <- LoadReference(path = "https://seurat.nygenome.org/azimuth/references/v1.0.0/human_pbmc")


## ----reading data---------------------------------------------------------------------------------
#query.data <- Read10X(data.dir = prefix_data)

#query.data <- readRDS(paste0(prefix_data, project.name, "_SingleCellExperiment.RDS"))
#query.data <- GetAssayData(object = query.data,
#                           slot = "data")
#query <- CreateSeuratObject(counts = query.data,
#                            project = project.name,
#                            min.cells = 3,
#                            min.features = 200)
#query


#barcodes <- read_tsv("./data/GSE132338/GSE132338_matrix_colData.txt.gz")
#barcodes <- barcodes %>%
#  select(donor)
#write_tsv(barcodes,
#          file = paste0("./data/GSE132338/barcodes.tsv"),
#          col_names = FALSE)


#features <- read_table(file = "./data/GSE132338/GSE132338_matrix_rowData.mtx.gz")
#write_tsv(features,
#          file = paste0("./data/GSE132338/features.tsv"),
#          col_names = FALSE)


#query.data <- Read10X(data.dir = prefix_data)


#query_raw <- read_table(file = "./data/GSE132338/GSE132338_matrix_counts.mtx.gz")


#query_raw <- read.table(file = "./data/GSE132338/GSE132338_matrix_counts.mtx.gz",
#                       as.is = TRUE)


query_raw <- readMM("./data/GSE132338/GSE132338_matrix_counts.mtx.gz")
barcodes <- read_tsv("./data/GSE132338/GSE132338_matrix_colData.txt.gz")
barcodes <- barcodes %>%
  dplyr::select(donor)
colnames(barcodes) <- NULL

features <- read_table(file = "./data/GSE132338/GSE132338_matrix_rowData.mtx.gz")
colnames(features) <- NULL

query_raw@Dimnames <- c(features, barcodes)

query <- CreateSeuratObject(counts = query_raw,
                            project = project.name,
                            min.cells = 3,
                            min.features = 200)
query

## ----loading batch_file and metadata--------------------------------------------------------------
query.colData <- read.table(paste0(prefix_data, project.name, "_matrix_colData.txt.gz"),
                        sep = "\t")
batch <- query.colData %>%
  dplyr::select(Cell = Barcode,
                batch = library.id)

batch_select <- query@assays[["RNA"]]@data@Dimnames[[2]] %>%
  tibble()
colnames(batch_select) <- c("Cell")
batch <- right_join(batch,
                    batch_select,
                    by = "Cell")

#metadata <- query.colData %>%
#  select(batch = library.id,
#         disease = classification) %>%
#  group_by(batch) %>%
#  arrange(batch) %>%
#  distinct(batch,
#           .keep_all = TRUE)
#write.csv(metadata,
#          file = paste0(prefix_data, project.name, "_metadata.csv"),
#          row.names = FALSE,
#          quote = FALSE)

#metadata_2 <- query.colData %>%
#  tibble() %>%
#  select(batch = library.id,
#         disease = classification,
#         celltype) %>%
#  group_by(batch,
#           celltype) %>%
#  dplyr::count()
#write.csv(metadata_2,
#          file = paste0(prefix_data, project.name, "_metadata.csv"),
#          row.names = FALSE,
#          quote = FALSE)

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

