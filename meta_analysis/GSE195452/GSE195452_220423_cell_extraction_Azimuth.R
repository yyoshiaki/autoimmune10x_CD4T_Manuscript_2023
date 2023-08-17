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
project.name <- c("GSE195452")
prefix_data <- paste0("./data/", project.name, "/")
prefix <- paste0("./output/", project.name, "/", project.name)


## ----Load_the_reference---------------------------------------------------------------------------
reference <- LoadReference(path = "https://seurat.nygenome.org/azimuth/references/v1.0.0/human_pbmc")


## ----reading data---------------------------------------------------------------------------------
df.1 <- read.delim('data/GSE195452/GSE195452-GPL18573_series_matrix.txt.gz', skip = 31, sep = "\t")
df.1 <- t(df.1[1:11,])
df.2 <- read.delim('data/GSE195452/GSE195452-GPL24676_series_matrix.txt.gz', skip = 31, sep = "\t")
df.2 <- t(df.2[1:11,])
df <- rbind(df.1, df.2)
df <- df[df[,9] == 'tissue: Blood',]
df <- as.data.frame(df)
df['PID'] <- df[,11] %>% str_remove("patient id: ")
df <- df[c(1,'PID')]
colnames(df) <- c('gsm', 'PID')
df <- as_tibble(df, rownames = "expid")

head(df)

df.metadata <- read_excel('data/GSE195452/1-s2.0-S0092867422003142-mmc1.xlsx', 
                     skip = 2) %>% 
  left_join(df, by = "PID") %>%
  drop_na(gsm) %>%
  rename(sample = PID, age = Age, sex = Gender, disease_duration = Duration, 
         disease = Disease,
         condition = Group, immunosuppressant = Therapy_DMARDs) %>%
  mutate(race = NA, immunosuppressant_duration = NA, batch = sample) %>%
  mutate(sex = fct_recode(sex, "Female" = "F", "Male" = "M"),
         disease = fct_recode(disease, "HC" = "Control", "SSc" = "SSC"))

write_csv(df.metadata %>% select(batch, sample, disease, disease_duration,
                                 condition, age, race, sex, immunosuppressant,
                                 immunosuppressant_duration), 
          paste0(prefix_data, project.name, "_metadata.csv"))

head(df.metadata)

list.q <- list()
for (r in rownames(df.metadata)) {
  gsm <- as.character(df.metadata[r,'gsm'])
  expid <- as.character(df.metadata[r,'expid'])
  s <- as.character(df.metadata[r,'sample'])
  print(s)
  mat <- read.table(paste0('data/',project.name,'/',gsm,'_',expid,'.txt.gz'))
  met <- as.data.frame(t(mat[c(),]))
  met['sample'] = s
  met['gsm'] = gsm
  met['expid'] = expid
  q <- CreateSeuratObject(mat, meta.data = met, project = s)
  list.q <- append(list.q, q)
}

query <- merge(list.q[[1]], y=list.q[2:length(list.q)],
               project = project.name)

# query@meta.data['sample'] <- query@meta.data['orig.ident']
d <- left_join(query@meta.data, as.data.frame(df.metadata))
rownames(d) <- rownames(query@meta.data)
query@meta.data <- d


## ----cell extraction-----------------------------------------------------------------------------------

extract_cells(query, reference, prefix)


## ----save.image-----------------------------------------------------------------------------------
prefix = paste0("./output/", project.name, "/")
save.image(file = paste0(prefix, "220423_cell_extraction_Azimuth.RData"))


## ----sessionInfo----------------------------------------------------------------------------------
sessionInfo()

