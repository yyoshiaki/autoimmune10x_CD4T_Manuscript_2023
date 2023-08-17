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
project.name <- c("GSE135779")
prefix_data <- paste0("./data/", project.name, "/")
prefix <- paste0("./output/", project.name, "/", project.name)


## ----Load_the_reference---------------------------------------------------------------------------
reference <- LoadReference(path = "https://seurat.nygenome.org/azimuth/references/v1.0.0/human_pbmc")


## ----reading data---------------------------------------------------------------------------------
df.samples <- read_csv('./data/GSE135779/GSE135779_samples.csv')
l.files <- list.files(paste0("./data/", project.name))
features <- read_table(file = "./data/GSE135779/GSE135779_genes.tsv.gz", col_names = FALSE)$X2

df.metadata <- read_csv("./data/GSE135779/41590_220_743_MOESM3_ESM_sheet1.csv")
df.metadata <- df.metadata %>% mutate(batch = Names, sample = Names, disease = Groups, 
                       disease_duration = NA, condition = paste0(Groups, ';SLEDAI:', SLEDAI),
                       age = Age, race = Race, sex = Gender, 
                       immunosuppressant = MMF, immunosuppressant_duration = NA) %>%
  mutate(disease = as.factor(disease) %>% fct_recode("HC"="aHD", "HC"="cHD", 
                                                     "SLE"="aSLE", "SLE"="cSLE"),
         sex = as.factor(sex) %>% fct_recode("Female"="F", "Male"="M")) %>%
  select(batch, sample, disease, disease_duration, condition, age, race, sex, 
           immunosuppressant, immunosuppressant_duration)
write_csv(df.metadata, paste0(prefix_data, project.name, "_metadata.csv")) 

for (r in rownames(df.samples)) {
  gsm <- as.character(df.samples[r,'GSM'])
  s <- as.character(df.samples[r,'sample'])
  prefix = paste0("./output/", project.name, "/", s)
  
  l.f <- l.files[startsWith(l.files, gsm)]
  f.b <- l.f[endsWith(l.f, "_barcodes.tsv.gz")]
  f.m <- l.f[endsWith(l.f, "_matrix.mtx.gz")]
  
  mat <- readMM(paste0("./data/", project.name,"/", f.m))
  coldata <- read.table(paste0("./data/", project.name,"/", f.b), sep="\t")
  colnames(mat) = coldata$V1
  rownames(mat) = features
  
  query <- CreateSeuratObject(counts = mat,
                              row.names = features,
                              project = project.name,
                              min.cells = 3,
                              min.features = 200)
  query@meta.data['sample'] <- s
  d <- merge(query@meta.data, df.metadata)
  rownames(d) <- row.names(query@meta.data)
  query@meta.data <- d
  extract_cells(query, reference, prefix)
}

## ----save.image-----------------------------------------------------------------------------------
prefix = paste0("./output/", project.name, "/")
save.image(file = paste0(prefix, "220423_cell_extraction_Azimuth.RData"))


## ----sessionInfo----------------------------------------------------------------------------------
sessionInfo()

