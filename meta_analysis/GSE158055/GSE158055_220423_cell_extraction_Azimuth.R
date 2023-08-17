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
project.name <- c("GSE158055")
prefix_data <- paste0("./data/", project.name, "/")
prefix <- paste0("./output/", project.name, "/", project.name)

## ----Load_the_reference---------------------------------------------------------------------------
reference <- LoadReference(path = "https://seurat.nygenome.org/azimuth/references/v1.0.0/human_pbmc")

## ----reading data---------------------------------------------------------------------------------
## ----reading data---------------------------------------------------------------------------------
barcode.path <- paste0(prefix_data, "barcodes.tsv.gz")
features.path <- paste0(prefix_data, "features.tsv.gz")
# matrix.path <- paste0(prefix_data, "matrix.mtx.gz")
# mat <- readMM(file = matrix.path)

### https://github.com/satijalab/seurat/issues/4030#issuecomment-776565733 (swapped $1 and $2)
# sed -n '4,1159554303p' GSE158055_covid19_counts.mtx |awk '{print $1"\t"$2"\t"$3}' > part1 &
# sed -n '1159554304,2319108602p' GSE158055_covid19_counts.mtx |awk '{print $1"\t"$2"\t"$3}' > part2
# vim part1
# vim part2

mat1 <- readMM(file = "./data/GSE158055/part1")
mat2 <- readMM(file = "./data/GSE158055/part2")

feature.names = read.delim(features.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)

colnames(mat1) = barcode.names$V1
rownames(mat1) = feature.names$V1
colnames(mat2) = barcode.names$V1
rownames(mat2) = feature.names$V1

q1 <- CreateSeuratObject(counts = mat1,
                         project = project.name,
                         assay = "RNA",
                         min.cells = 3,
                         min.features = 200)

q2 <- CreateSeuratObject(counts = mat2,
                         project = project.name,
                         assay = "RNA", 
                         min.cells = 3,
                         min.features = 200)

intersect(Cells(q1), Cells(q2))
## list the cells (=colnames) that you want to remove:
toRemove <- c("GATGAAAAGTCCTCCT-151")
## filter them out:
q1_filtered <- q1[,!colnames(q1) %in% toRemove]
intersect(Cells(q1_filtered), Cells(q2))

#q <- merge(q1_filtered,
#           y = q2,
#           project = project.name)

cells <- read_csv(paste0(prefix_data, "GSE158055_cell_annotation.csv.gz"))
cells.1 <- tibble(cells$sampleID)
rownames(cells.1) <- cells$cellName
colnames(cells.1) <- c("sample_name")
q1_filtered <- AddMetaData(q1_filtered, cells.1)
q2 <- AddMetaData(q2, cells.1)

meta_s <- read_csv(paste0(prefix_data, "GSE158055_metadata_sample.csv"))
q1_filtered <- subset(q1_filtered, subset = sample_name %in% meta_s$sample_name)
q2 <- subset(q2, subset = sample_name %in% meta_s$sample_name)

cells.2 <-  tibble(cells$cellName, cells$sampleID)
colnames(cells.2) <- c("cellName", "sample_name")
cells.2 <- subset(cells.2, subset = sample_name %in% meta_s$sample_name)

meta <- read_csv(paste0(prefix_data, "GSE158055_metadata.csv"))
cells.2 <- cells.2 %>%
  group_by(sample_name) %>%
  left_join(meta_s, by = c("sample_name"))
cells.2 <- cells.2 %>%
  group_by(sample) %>%
  left_join(meta, by = c("sample"))

cells.q1_filtered <- subset(cells.2,
                            subset = sample_name %in% unique(q1_filtered@meta.data[["sample_name"]]))
cells.q2 <- subset(cells.2,
                   subset = sample_name %in% unique(q2@meta.data[["sample_name"]]))

rownames(cells.q1_filtered) <- cells.q1_filtered$cellName
rownames(cells.q2) <- cells.q2$cellName

q1_filtered <- AddMetaData(q1_filtered, cells.q1_filtered)
q2 <- AddMetaData(q2, cells.q2)


## ----process--------------------------------------------------------------------------------------
obj.list.q1_filtered <- SplitObject(q1_filtered, split.by = "sample")
ss_1 <- 0
for (i in unique(q1_filtered@meta.data[["sample"]])) {
  print(i)
  if (obj.list.q1_filtered[[i]]@assays[["RNA"]]@counts@Dim[2] >= 500) {
    prefix <- paste0("./output/", project.name, "/", i)
    extract_cells(obj.list.q1_filtered[[i]], reference, prefix)
    if (ss_1==0) {
      ii <- tibble(i)
      colnames(ii) <- "sample"
      ss_1 <- ii
    } else {
      ii <- tibble(i)
      colnames(ii) <- "sample"
      ss_1 <- ss_1 %>% bind_rows(ii)
    }
  } else {
    print(paste0(i, " was removed from processing because cell number is <500"))
  }
}


obj.list.q2 <- SplitObject(q2, split.by = "sample")
ss_2 <- 0
for (i in unique(q2@meta.data[["sample"]])) {
  print(i)
  if (obj.list.q2[[i]]@assays[["RNA"]]@counts@Dim[2] >= 500) {
    prefix <- paste0("./output/", project.name, "/", i)
    extract_cells(obj.list.q2[[i]], reference, prefix)
    if (ss_2==0) {
      ii <- tibble(i)
      colnames(ii) <- "sample"
      ss_2 <- ii
    } else {
      ii <- tibble(i)
      colnames(ii) <- "sample"
      ss_2 <- ss_2 %>% bind_rows(ii)
    }
  } else {
    print(paste0(i, " was removed from processing because cell number is <500"))
  }
}

ss <- ss_1 %>%
  bind_rows(ss_2)
write_csv(ss, file = paste0("./output/", project.name, "/", project.name, "_processed_sample_cell_extraction.csv"))


## ----save.image-----------------------------------------------------------------------------------
#prefix = paste0("./output/", project.name, "/")
#save.image(file = paste0(prefix, "220423_cell_extraction_Azimuth.RData"))


## ----sessionInfo----------------------------------------------------------------------------------
sessionInfo()

