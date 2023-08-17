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
project.name <- c("GSE168732")
prefix_data <- paste0("./data/", project.name, "/")
prefix <- paste0("./output/", project.name, "/", project.name)

## ----Load_the_reference---------------------------------------------------------------------------
reference <- LoadReference(path = "https://seurat.nygenome.org/azimuth/references/v1.0.0/human_pbmc")

## ----reading data---------------------------------------------------------------------------------
df.metadata <- read_csv(paste0(prefix_data, project.name, "_metadata.csv"))

holder <- list()
for (s in df.metadata$sample){
  mat <- readMM(paste0("./data/",project.name,"/", s ,"_matrix.mtx.gz"))
  barcodes <- read_tsv(paste0("./data/",project.name,"/", s ,"_barcodes.tsv.gz"),
                       col_names = FALSE)$X1
  features <- read_tsv(file = paste0("./data/",project.name,"/", s ,"_features.tsv.gz"),
                       col_names = FALSE)$X2
  colnames(mat) = barcodes
  rownames(mat) = features
  
  q <- CreateSeuratObject(counts = mat,
                              project = s,
                              min.cells = 3,
                              min.features = 200)
  holder <- c(holder, q)
}

query <- merge(holder[[1]], y=holder[2:length(holder)],
               add.cell.ids = df.metadata$sample,
               project = project.name)

d <- merge(query@meta.data, df.metadata, by.x = "orig.ident", by.y = "sample")
row.names(d) <- row.names(query@meta.data)
query@meta.data <- d
query@meta.data['sample'] <- query@meta.data['orig.ident']

## ----process--------------------------------------------------------------------------------------
extract_cells(query, reference, prefix)


## ----save.image-----------------------------------------------------------------------------------
prefix = paste0("./output/", project.name, "/")
save.image(file = paste0(prefix, "220423_cell_extraction_Azimuth.RData"))


## ----sessionInfo----------------------------------------------------------------------------------
sessionInfo()
