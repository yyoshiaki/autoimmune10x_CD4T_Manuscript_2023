#install.packages("anndata")


## ----library and source---------------------------------------------------------------------------
library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(anndata)
library(Matrix)

setwd('/home/rstudio/autoimmune_10x')



## ----parameter setting----------------------------------------------------------------------------
project.name <- c("GSE174188")
prefix_data <- paste0("./data/", project.name, "/")
prefix <- paste0("./output/", project.name, "/", project.name)




## ----reading data---------------------------------------------------------------------------------
#Convert(paste0(prefix_data, "GSE174188_CLUES1_adjusted.h5ad"), dest = "h5seurat", overwrite = TRUE)
cells <- read_h5ad(paste0(prefix_data, "GSE174188_CLUES1_adjusted.h5ad"))
cells <- as(object = cells, Class = "dgCMatrix")
cells <- CreateSeuratObject(counts = t(cells$X),
                            row.names = cells$var,
                            meta.data = cells$obs)
saveRDS(cells,
        file = paste0(prefix_data, "GSE174188_CLUES1_adjusted_seuratobj.rds"))



## ----sessionInfo----------------------------------------------------------------------------------
sessionInfo()

