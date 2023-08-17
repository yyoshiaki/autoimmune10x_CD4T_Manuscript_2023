## ----Ensure, eval=FALSE, include=FALSE------------------------------------------------------------------------
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


## ----library--------------------------------------------------------------------------------------------------
library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(Azimuth)
library(patchwork)
library(tidyverse)
library(sctransform)

setwd('/home/rstudio/autoimmune_10x')
source('./script/functions.R')

## ----parameter setting----------------------------------------------------------------------------------------
project.name <- "E-MTAB-10026"
prefix_data <- paste0("./data/", project.name, "/")
prefix <- paste0("./output/", project.name, "/", project.name)

## ----Load_the_reference---------------------------------------------------------------------------------------
reference <- LoadReference(path = "https://seurat.nygenome.org/azimuth/references/v1.0.0/human_pbmc")


## ----reading data---------------------------------------------------------------------------------------------
# Convert("./data/E-MTAB-10026/haniffa21.processed.h5ad", dest = "h5seurat", overwrite = TRUE)
query <- LoadH5Seurat("./data/E-MTAB-10026/haniffa21.processed.h5seurat")

query <- subset(query, subset = Status %in% c("Covid", "Healthy"))
query@meta.data$Status <- factor(as.character(query@meta.data$Status))
levels(query@meta.data$Status) <- list(Covid  = "Covid", HC = "Healthy")
query$sample_id <- factor(query$sample_id)
levels(query@meta.data$Age_interval) <- c(25, 35, 45, 55, 65, 75, 85, 95)
query@meta.data$Age_interval <- query@meta.data$Age_interval %>% as.character() %>% as.numeric()
query@meta.data <- query@meta.data %>% 
  mutate(sample = patient_id,
         disease = Status,
         disease_duration = NA,
         condition = Status_on_day_collection_summary,
         age = Age_interval,
         race = NA,
         sex = Sex,
         immunosuppressant = NA,
         immunosuppressant_duration = NA,
         batch = patient_id)
query@meta.data['project'] <- project.name
# query@meta.data['Cell'] <- rownames(query@meta.data)

query@meta.data <- query@meta.data %>%
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

cols <- c("batch", "sample", "disease", "disease_duration", "condition", "age", "race", "sex", 
"immunosuppressant","immunosuppressant_duration")
df.metadata <- query@meta.data[cols] %>%
  distinct()
write_csv(df.metadata, paste0(prefix_data, project.name, "_metadata.csv")) 

## ----process--------------------------------------------------------------------------------------------------
for (sample in levels(query@meta.data$sample)) {
  print(sample)
  prefix = paste0("./output/", project.name, "/", sample)
  
  s = sample
  if (sum(query@meta.data['sample'] == s) > 0) {
    q <- subset(query, subset = sample == as.character(s))
    q <- CreateSeuratObject(counts = q[["raw"]],
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


## ----save.image-----------------------------------------------------------------------------------------------
prefix = paste0("./output/", project.name, "/")
save.image(file = paste0(prefix, "220423_CD4Tcell_extraction_Azimuth.RData"))


## ----sessionInfo----------------------------------------------------------------------------------------------
sessionInfo()

