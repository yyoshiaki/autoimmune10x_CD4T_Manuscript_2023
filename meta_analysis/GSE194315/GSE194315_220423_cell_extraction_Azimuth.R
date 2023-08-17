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
project.name <- c("GSE194315")
prefix_data <- paste0("./data/", project.name, "/")
prefix <- paste0("./output/", project.name, "/", project.name)

## ----Load_the_reference---------------------------------------------------------------------------
reference <- LoadReference(path = "https://seurat.nygenome.org/azimuth/references/v1.0.0/human_pbmc")

## ----reading data---------------------------------------------------------------------------------
df.metadata <- read_excel('data/GSE194315/Table_1_Combined Single Cell Transcriptome and Surface Epitope Profiling Identifies Potential Biomarkers of Psoriatic Arthritis and Facilitates Diagno.xlsx', 
           sheet=2)

df.metadata <- df.metadata %>% mutate(sample = SubjectID, disease = Status, 
                        disease_duration = NA, condition = Status,
                        age = Age, race = Ethnicity,
                        sex = Gender, immunosuppressant = `Systemic Medications`,
                        immunosuppressant_duration = NA, batch = SubjectID) %>%
  mutate(disease = as.factor(disease)%>%
           fct_recode("HC" = "Healthy", "Psoriasis" = "PSO", 
                      "Psoriasis" = "PSA", "Psoriasis" = "PSX"),
         sex = as.factor(sex)%>%
           fct_recode("Female" = "F", "Male" = "M"),
         race = as.factor(race)%>%
           fct_recode("European" = "White")) %>%
  select(sample,
           disease,
           disease_duration,
           condition,
           age,
           race,
           sex,
           immunosuppressant,
           immunosuppressant_duration,
           batch)

write_csv(df.metadata, paste0(prefix_data, project.name, "_metadata.csv")) 

files <- list.files(paste0("./data/",project.name), pattern = "matrix.mtx.gz") %>%
  str_replace(".matrix.mtx.gz", "")

holder <- list()
for (f in files){
  mat <- readMM(paste0("./data/",project.name,"/", f ,".matrix.mtx.gz"))
  barcodes <- read_tsv(paste0("./data/",project.name,"/", f ,".barcodes.tsv.gz"),
                       col_names = FALSE)$X1
  features <- read_tsv(file = paste0("./data/",project.name,"/", f ,".features.tsv.gz"),
                       col_names = FALSE)$X2
  colnames(mat) = barcodes
  rownames(mat) = features
  
  q <- CreateSeuratObject(counts = mat,
                          project = f,
                          min.cells = 3,
                          min.features = 200)
  holder <- c(holder, q)
}

query <- merge(holder[[1]], y=holder[2:length(holder)],
               add.cell.ids = files,
               project = project.name)

# df.cellmeta.1 <- read_tsv('data/GSE194315/GSE194315_CellMetadata-PSA_TotalCiteseq_20220103.tsv.gz')
# df.cellmeta.2 <- read_tsv('data/GSE194315/GSE194315_CellMetadata-AS_TotalCiteseq_20220103.tsv.gz')
# df.cellmeta <- rbind(df.cellmeta.1, df.cellmeta.2)
df.cellmeta <- read_tsv('data/GSE194315/GSE194315_CellMetadata-PSA_TotalCiteseq_20220103.tsv.gz')

query@meta.data['CellName'] <- rownames(query@meta.data) %>% str_replace("-1$", "")
d <- query@meta.data
d <- left_join(d, as.data.frame(df.cellmeta))
d <- left_join(d, as.data.frame(df.metadata), by=c("Subject" = "sample"))
rownames(d) <- rownames(query@meta.data)
query@meta.data <- d

query <- subset(x = query, disease %in% c('HC', 'Psoriasis'))
query <- subset(x = query, subset = DemuxletDropletType == "SNG")
query@meta.data['sample'] <- query@meta.data['Subject']

## ----process--------------------------------------------------------------------------------------
extract_cells(query, reference, prefix)


## ----save.image-----------------------------------------------------------------------------------
prefix = paste0("./output/", project.name, "/")
save.image(file = paste0(prefix, "220423_cell_extraction_Azimuth.RData"))


## ----sessionInfo----------------------------------------------------------------------------------
sessionInfo()
