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
project.name <- c("1M_PBMC_T1D_Parse")
prefix_data <- paste0("./data/", project.name, "/")
prefix <- paste0("./output/", project.name, "/", project.name)

## ----Load_the_reference---------------------------------------------------------------------------
reference <- LoadReference(path = "https://seurat.nygenome.org/azimuth/references/v1.0.0/human_pbmc")

## ----reading data---------------------------------------------------------------------------------
barcode.path <- paste0(prefix_data, "cell_metadata_1M_PBMC.csv")
features.path <- paste0(prefix_data, "all_genes_1M_PBMC.csv")
matrix.path <- paste0(prefix_data, "DGE_1M_PBMC.mtx")
mat <- readMM(file = matrix.path)
feature.names = read.csv(features.path,
                           header = TRUE,
                           stringsAsFactors = FALSE)
barcode.names = read.csv(barcode.path,
                           header = TRUE,
                           stringsAsFactors = FALSE)

mat <- t(mat)
colnames(mat) <- barcode.names$bc_wells
rownames(mat) = feature.names$gene_name

mat <- mat[rownames(mat) != "",]

query <- CreateSeuratObject(counts = mat,
                            project = project.name,
                            assay = "RNA",
                            min.cells = 3,
                            min.features = 200)

metadata <- query@assays[["RNA"]]@counts@Dimnames[[2]] %>%
  tibble()
colnames(metadata) <- c("bc_wells")

meta <- read_csv(paste0(prefix_data, project.name, "_metadata.csv"))

metadata <- metadata %>% left_join(barcode.names, by = "bc_wells") %>%
  left_join(meta, by="sample")

row.names(metadata) <- metadata$bc_wells

query <- AddMetaData(query, metadata)

# filter omitted samples
query <- query[,!is.na(query@meta.data$disease)]

## ----process--------------------------------------------------------------------------------------------------
for (sample in unique(query@meta.data$sample)) {
  print(sample)
  prefix = paste0("./output/", project.name, "/", sample)
  
  s = sample
  if (sum(query@meta.data['sample'] == s) > 0) {
    q <- subset(query, subset = sample == as.character(s))
    q <- CreateSeuratObject(counts = q@assays[["RNA"]]@counts,
                            project = project.name,
                            assay = "RNA",
                            meta.data = q@meta.data,
                            min.cells = 3,
                            min.features = 200) 
    extract_cells(q, reference, prefix)
  } else {
    print('No cell included in the sample. skipped...')
  }
}



## ----save.image-----------------------------------------------------------------------------------
prefix = paste0("./output/", project.name, "/")
save.image(file = paste0(prefix, "220423_cell_extraction_Azimuth.RData"))


## ----sessionInfo----------------------------------------------------------------------------------
sessionInfo()

