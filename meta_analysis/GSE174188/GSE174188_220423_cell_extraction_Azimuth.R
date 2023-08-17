## ----library and source---------------------------------------------------------------------------
library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(Azimuth)
library(patchwork)
library(tidyverse)
library(sctransform)
library(Matrix)

setwd('/home/rstudio/autoimmune_10x')
source('./script/functions.R')


## ----parameter setting----------------------------------------------------------------------------
project.name <- c("GSE174188")
prefix_data <- paste0("./data/", project.name, "/")
prefix <- paste0("./output/", project.name, "/", project.name)


## ----Load_the_reference---------------------------------------------------------------------------
reference <- LoadReference(path = "https://seurat.nygenome.org/azimuth/references/v1.0.0/human_pbmc")


## ----reading data---------------------------------------------------------------------------------
Convert(paste0(prefix_data, "GSE174188_CLUES1_adjusted.raw.h5ad"),
        dest = "h5seurat",
        overwrite = TRUE)

so <- LoadH5Seurat(paste0(prefix_data, "GSE174188_CLUES1_adjusted.raw.h5seurat"),
                   meta.data = FALSE)
metadata <- read.csv(paste0(prefix_data, "obs.csv"),
                     row.names = 1)
metadata <- metadata %>%
  mutate(batch = ind_cov,
         sample = ind_cov,
         disease = gsub(SLE_status, pattern = "Healthy", replacement = "HC"),
         disease_duration = NA,
         condition = Status,
         age = Age,
         race = pop_cov,
         sex = str_replace_all(Sex, pattern = c("F"="f", "M"="m")),
         immunosuppressant = NA,
         immunosuppressant_duration = NA) %>%
  relocate(sample,
           disease,
           disease_duration,
           condition,
           age,
           race,
           sex,
           immunosuppressant,
           immunosuppressant_duration,
           batch)
  
so@meta.data <- metadata

options(repr.plot.width = 10, repr.plot.height = 8)
FeaturePlot(object = so,
            features = c("CD3E", 'FOXP3', 'MX1'),
            reduction = "umap",
            ncol = 2,
            raster = TRUE)
ggsave(paste0(prefix, "_featureplot_raw_h5seurat.pdf"),
       width = 10,
       height = 10)

options(repr.plot.height = 4, repr.plot.width = 10)
DimPlot(so,
        reduction = "umap",
        group.by = "SLE_status",
        shuffle = TRUE) +
  labs(title = "SLE_status")
ggsave(paste0(prefix, "_dimplot_raw_h5seurat.pdf"),
       width = 5,
       height = 5)



obj <- CreateSeuratObject(counts = so@assays$RNA@counts,
                          project = project.name,
                          assay = "RNA",
                          min.cells = 3,
                          min.features = 200)

obj.meta.data <- so@meta.data
obj <- AddMetaData(obj,
                   metadata = obj.meta.data)

#m <- metadata %>%
#  group_by(sample)  %>%
#  distinct(ind_cov_batch_cov, .keep_all = TRUE) %>%
#  arrange(sample)
#m <- ungroup(m)
#m %>%
#  count(condition)

print(paste0("---Before the removal of Treated sample---"))
m <- obj@meta.data
m <- m %>%
  group_by(sample)  %>%
  distinct(ind_cov_batch_cov, .keep_all = TRUE) %>%
  arrange(sample)
m <- ungroup(m)
count(m, condition)
sum(count(m, condition)[2])
nrow(count(m, sample))
print(paste0("---Before the removal of Treated sample---"))


print(paste0("---After the removal of Treated sample---"))
obj <- subset(obj,
              subset = condition!=c("Treated"))
m <- obj@meta.data
m <- m %>%
  group_by(sample)  %>%
  distinct(ind_cov_batch_cov, .keep_all = TRUE) %>%
  arrange(sample)
m <- ungroup(m)
count(m, condition)
sum(count(m, condition)[2])
nrow(count(m, sample))
print(paste0("---After the removal of Treated sample---"))



m <- m %>%
  group_by(sample) %>%
  distinct(sample)
nrow(m)
write_csv(m, file = paste0(prefix, "_sample.csv"))

obj.list <- SplitObject(obj, split.by = "sample")
for (i in m$sample) {
  print(i)
  prefix <- paste0("./output/", project.name, "/", i)
  extract_cells(obj.list[[i]], reference, prefix)
}


## ----save.image-----------------------------------------------------------------------------------
prefix = paste0("./output/", project.name, "/")
save.image(file = paste0(prefix, "220423_cell_extraction_Azimuth.RData"))


## ----sessionInfo----------------------------------------------------------------------------------
sessionInfo()

