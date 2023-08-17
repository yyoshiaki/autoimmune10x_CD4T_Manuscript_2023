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

setwd('/home/rstudio/autoimmune_10x')
source('./script/functions.R')


## ----parameter setting----------------------------------------------------------------------------
project.name <- c("GSE144744")
prefix_data <- paste0("./data/", project.name, "/")
prefix <- paste0("./output/", project.name, "/", project.name)

## ----Load_the_reference---------------------------------------------------------------------------
reference <- LoadReference(path = "https://seurat.nygenome.org/azimuth/references/v1.0.0/human_pbmc")

## ----reading data---------------------------------------------------------------------------------
#metadata <- read_csv(paste0(prefix_data, "GSE144744_metadata_per_sample.csv.gz"))
metadata <- read_csv(paste0(prefix_data, "GSE144744_metadata_per_cell.csv.gz"))
meta <- read_csv(paste0(prefix_data, "GSE144744_metadata.csv"))
metadata <- metadata %>%
  select(cell_names, sample_10X) %>%
  group_by(sample_10X) %>%
  left_join(meta, by = "sample_10X")
row.names(metadata) <- metadata$cell_names


matrix_dir = "./data/GSE144744/RNA_counts/"
barcode.path <- paste0(matrix_dir, "barcodes.tsv")
features.path <- paste0(matrix_dir, "genes.tsv")
matrix.path <- paste0(matrix_dir, "matrix.mtx")
mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$V1


query <- CreateSeuratObject(counts = mat,
                            project = project.name,
                            assay = "RNA",
                            min.cells = 3,
                            min.features = 200)

query <- AddMetaData(query, metadata)

query$orig.ident = factor(paste(project.name))
#Idents(query) <- "orig.ident"


## ----process--------------------------------------------------------------------------------------
sample <- meta %>%
  select(sample) %>%
  group_by(sample) %>%
  distinct()


for (s in sample$sample) {
  print(s)
  prefix = paste0("./output/", project.name, "/", s)
  
  q <- query %>%
    subset(sample == s)
  
  extract_cells(q, reference, prefix)
  
}


## ----save.image-----------------------------------------------------------------------------------
prefix = paste0("./output/", project.name, "/")
save.image(file = paste0(prefix, "220423_cell_extraction_Azimuth.RData"))


## ----sessionInfo----------------------------------------------------------------------------------
sessionInfo()
