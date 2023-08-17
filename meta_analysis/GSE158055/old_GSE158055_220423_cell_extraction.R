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
project.name <- c("GSE158055")
prefix_data <- paste0("./data/", project.name, "/")
prefix <- paste0("./output/", project.name, "/", project.name)

## ----Load_the_reference---------------------------------------------------------------------------
reference <- LoadReference(path = "https://seurat.nygenome.org/azimuth/references/v1.0.0/human_pbmc")

## ----reading data---------------------------------------------------------------------------------
barcode.path <- paste0(prefix_data, "barcodes.tsv.gz")
features.path <- paste0(prefix_data, "features.tsv.gz")
# matrix.path <- paste0(prefix_data, "matrix.mtx.gz")
# mat <- readMM(file = matrix.path)

### https://github.com/satijalab/seurat/issues/4030#issuecomment-776565733 (swapped $1 and $2)
# sed -n '4,1159554303p' GSE158055_covid19_counts.mtx |awk '{print $1"\t"$2"\t"$3}' > part1 &
# sed -n '1159554304,2319108602p' GSE158055_covid19_counts.mtx |awk '{print $1"\t"$2"\t"$3}' > part2
# vim part1
# vim part2

mat1 <- readMM(file = "./data/GSE158055/part1")
mat2 <- readMM(file = "./data/GSE158055/part2")

feature.names = read.delim(features.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)

colnames(mat1) = barcode.names$V1
rownames(mat1) = feature.names$V1
colnames(mat2) = barcode.names$V1
rownames(mat2) = feature.names$V1

q1 <- CreateSeuratObject(counts = mat1,
                            project = project.name,
                            assay = "RNA",
                            min.cells = 3,
                            min.features = 200)

q2 <- CreateSeuratObject(counts = mat2,
                         project = project.name,
                         assay = "RNA", 
                         min.cells = 3,
                         min.features = 200)

metadata <- query@assays[["RNA"]]@counts@Dimnames[[2]] %>%
  tibble()
colnames(metadata) <- c("labels")

meta_s <- read_csv(paste0(prefix_data, "GSE149689_metadata_sample.csv"))
meta_s$index <- as.character(meta_s$index)
meta <- read_csv(paste0(prefix_data, "GSE149689_metadata.csv"))
metadata <- metadata %>%
  mutate(index = str_split(string = metadata$labels, pattern = "-", n = 2, simplify = TRUE)[,2])
metadata <- metadata %>%
  group_by(index) %>%
  left_join(meta_s, by = c("index"))
metadata <- metadata %>%
  group_by(sample) %>%
  left_join(meta, by = c("sample"))
row.names(metadata) <- metadata$labels

query <- AddMetaData(query, metadata)


## ----process--------------------------------------------------------------------------------------
extract_cells(query, reference, prefix)


## ----save.image-----------------------------------------------------------------------------------
prefix = paste0("./output/", project.name, "/")
save.image(file = paste0(prefix, "220423_cell_extraction_Azimuth.RData"))


## ----sessionInfo----------------------------------------------------------------------------------
sessionInfo()

