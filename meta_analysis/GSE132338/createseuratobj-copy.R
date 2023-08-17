library(tidyverse)
library(Seurat)
library(Matrix)

setwd('~/autoimmune_10x/')

project.name <- c("GSE132338")

mat <- readMM("./data/GSE132338/GSE132338_matrix_counts.mtx.gz")
coldata <- read.table("./data/GSE132338/GSE132338_matrix_colData.txt.gz", sep="\t")
features <- read_table(file = "./data/GSE132338/GSE132338_matrix_rowData.mtx.gz")
colnames(mat) = rownames(coldata)
rownames(mat) = features$gene.symbol

query <- CreateSeuratObject(counts = mat,
                              row.names = features$gene.symbol,
                              meta.data = coldata,
                              project = project.name,
                              min.cells = 3,
                              min.features = 200)
