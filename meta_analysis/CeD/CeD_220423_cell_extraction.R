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
library(readxl)

source('./script/functions.R')


## ----parameter setting----------------------------------------------------------------------------
project.name <- c("CeD")
prefix_data <- paste0("./data/", project.name, "/")
prefix <- paste0("./output/", project.name, "/", project.name)

## ----Load_the_reference---------------------------------------------------------------------------
reference <- LoadReference(path = "https://seurat.nygenome.org/azimuth/references/v1.0.0/human_pbmc")

## ----reading data---------------------------------------------------------------------------------

mat <- read.csv('data/CeD/Supplementary_Table_8_scRNAseq_rawcounts.csv')
metadata <- read.csv('data/CeD/DataSheet_1_Single-Cell RNA Sequencing of Peripheral Blood Mononuclear Cells From Pediatric Coeliac Disease Patients Suggests Potential Pre-Seroconver.csv')
rownames(metadata) <- colnames(mat)

query <- CreateSeuratObject(counts = mat,
                        project = project.name,
                        meta.data = metadata,
                        min.cells = 3,
                        min.features = 200)

query@meta.data['cond_time'] = paste0(query@meta.data$Condition, ";", query@meta.data$Sampling_time)

query <- subset(query, cond_time != "CeD;T1")

query@meta.data <- query@meta.data %>% mutate(sample = paste0(Paring_code, ';', Condition),
                           disease = as.factor(Condition)%>%
                             fct_recode("HC" = "Controls", "Celiac" = "CeD"),
                           disease_duration = NA, 
                           condition = as.factor(Condition)%>%
                             fct_recode("HC" = "Controls", "Celiac" = "CeD"),
                           age = Age_at_sampling_months/12,
                           race = NA,
                           sex = as.factor(Sex)%>%
                             fct_recode("Female" = "F", "Male" = "M"),
                           immunosuppressant = 0, immunosuppressant_duration = NA, 
                           batch = Paring_code)

df.metadata <-  query@meta.data %>%
  select(sample,
         disease,
         disease_duration,
         condition,
         age,
         race,
         sex,
         immunosuppressant,
         immunosuppressant_duration,
         batch) %>%
  distinct() %>% 
  remove_rownames()

write_csv(df.metadata, paste0(prefix_data, project.name, "_metadata.csv")) 

## ----process--------------------------------------------------------------------------------------
extract_cells(query, reference, prefix)


## ----save.image-----------------------------------------------------------------------------------
prefix = paste0("./output/", project.name, "/")
save.image(file = paste0(prefix, "220423_cell_extraction_Azimuth.RData"))


## ----sessionInfo----------------------------------------------------------------------------------
sessionInfo()
