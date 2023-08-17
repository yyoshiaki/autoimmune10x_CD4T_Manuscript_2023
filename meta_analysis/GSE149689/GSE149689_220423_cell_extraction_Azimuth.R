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
project.name <- c("GSE149689")
prefix_data <- paste0("./data/", project.name, "/")
prefix <- paste0("./output/", project.name, "/", project.name)

## ----Load_the_reference---------------------------------------------------------------------------
reference <- LoadReference(path = "https://seurat.nygenome.org/azimuth/references/v1.0.0/human_pbmc")

## ----reading data---------------------------------------------------------------------------------
barcode.path <- paste0(prefix_data, "barcodes.tsv.gz")
features.path <- paste0(prefix_data, "features.tsv.gz")
matrix.path <- paste0(prefix_data, "matrix.mtx.gz")
mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$V2


query <- CreateSeuratObject(counts = mat,
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

