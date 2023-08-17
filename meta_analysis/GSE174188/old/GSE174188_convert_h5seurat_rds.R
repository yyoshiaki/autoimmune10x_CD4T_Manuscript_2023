## ----library and source---------------------------------------------------------------------------
library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(tidyverse)

setwd('/home/rstudio/autoimmune_10x')



## ----parameter setting----------------------------------------------------------------------------
project.name <- c("GSE174188")
prefix_data <- paste0("./data/", project.name, "/")
prefix <- paste0("./output/", project.name, "/", project.name)




## ----reading data---------------------------------------------------------------------------------
Convert(paste0(prefix_data, "GSE174188_CLUES1_adjusted.raw.h5ad"),
        dest = "h5seurat",
        overwrite = TRUE)

so <- LoadH5Seurat(paste0(prefix_data, "GSE174188_CLUES1_adjusted.raw.h5seurat"),
                   meta.data = FALSE)
metadata <- read.csv(paste0(prefix_data, "obs.csv"),
                     row.names = 1)

so@meta.data <- metadata

options(repr.plot.width = 10, repr.plot.height = 8)
FeaturePlot(object = so,
            features = c("CD3E", 'FOXP3', 'MX1'),
            reduction = "umap",
            ncol = 2,
            raster = TRUE)
ggsave(paste0(prefix_data, project.name, "_featureplot_raw_h5seurat.pdf"),
       width = 10,
       height = 10)

options(repr.plot.height = 4, repr.plot.width = 10)
DimPlot(so,
        reduction = "umap",
        group.by = "SLE_status",
        shuffle = TRUE) +
  labs(title = "SLE_status")
ggsave(paste0(prefix_data, project.name, "_dimplot_raw_h5seurat.pdf"),
       width = 5,
       height = 5)

saveRDS(so,
        file = paste0(prefix_data, "GSE174188_CLUES1_adjusted.raw.rds"))


## ----sessionInfo----------------------------------------------------------------------------------
sessionInfo()

