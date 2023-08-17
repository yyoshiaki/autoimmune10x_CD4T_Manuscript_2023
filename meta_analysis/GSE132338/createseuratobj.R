library(tidyverse)
library(Seurat)
library(Matrix)

setwd('/home/dtakeuchi/autoimmune_pbmc_10x')

# query_raw <- read_table(file = "./data/GSE132338/GSE132338_matrix_counts.mtx.gz")

query_raw <- readMM("./data/GSE132338/GSE132338_matrix_counts.mtx.gz")
coldata <- read.table("./data/GSE132338/GSE132338_matrix_colData.txt.gz", sep="\t")


features <- read_table(file = "./data/GSE132338/GSE132338_matrix_rowData.mtx.gz")

query <- CreateSeuratObject(counts = query_raw,
                              row.names = features,
                              meta.data = coldata,
                              project = project.name,
                              min.cells = 3,
                              min.features = 200)
