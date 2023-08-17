library(tidyverse)
library(ReactomePA)
library(clusterProfiler)

library(pheatmap)
library(pals)

setwd('/home/yyasumizu/media32TB/bioinformatics/autoimmune_10x/R')
df <- read_csv('../scanpy/220423_STR4.5.21_MS01.02.03.04_MG01.02.03_SL02.03.04_Tcell/NMF_corr_top100.csv')
df

n_top <- 50

list.genes <- list()
for (c in 0:11) {
  print(c)
  x <- df[1:n_top, paste0('gene_NMF', c)]  %>% unlist(use.names = FALSE)
  eg <- bitr(x, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  list.genes[[paste0('NMF', c)]] <-  eg$ENTREZID
}

ck <- compareCluster(geneCluster = list.genes, fun = enrichPathway)
# ck <- setReadable(ck, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
head(ck) 

dotplot(ck)
ggsave("../scanpy/220423_STR4.5.21_MS01.02.03.04_MG01.02.03_SL02.03.04_Tcell/graph/reactome.dotplot.pdf", width = 15, height = 12)

write.csv(ck, '../scanpy/220423_STR4.5.21_MS01.02.03.04_MG01.02.03_SL02.03.04_Tcell/NMF_reactome_top50.csv')

# row.labels <- c('NMF0 Cytotoxic', 'NMF1 Treg', 'NMF2 Th17', 'NMF3 Naiveness', 
#                 'NMF4 Act', 'NMF5 Th2', 'NMF6 Tfh', 'NMF7 IFN', 'NMF8 Cent. Mem.',
#                 'NMF9 Thymic Emi.', 'NMF10 Tissue', 'NMF11 Th1')
row.labels <- c('NMF0', 'NMF1', 'NMF2', 'NMF3', 
                'NMF4', 'NMF5', 'NMF6', 'NMF7', 'NMF8',
                'NMF9', 'NMF10', 'NMF11')
n_top <- 15
list.genes <- list()
for (c in 0:11) {
  print(c)
  x <- df[1:n_top, paste0('gene_NMF', c)]  %>% unlist(use.names = FALSE)
  list.genes[[paste0('NMF', c)]] <- x
}

df.corr <- read.csv('../scanpy/220438_NMFprojection/NMF_corr_all.csv', row.names = 1)
genes <- c(list.genes$NMF0, list.genes$NMF1, list.genes$NMF2, list.genes$NMF3, list.genes$NMF4,
           list.genes$NMF5, list.genes$NMF6, list.genes$NMF7, list.genes$NMF8, list.genes$NMF9,
           list.genes$NMF10, list.genes$NMF11)
d.plot <- df.corr[genes,]

p <- pheatmap(t(d.plot), cluster_rows = FALSE, cluster_cols = FALSE,
         gaps_row=seq(0,11), gaps_col = seq(0, 12*n_top-1, n_top),
         color =ocean.balance(100), breaks = seq(-1, 1, by=0.02), labels_col=genes)
p

pdf('../scanpy/220438_NMFprojection/graph/pheatmap.corr.pdf', width = 24, height = 3)
p
dev.off()

n_top <- 10
list.genes <- list()
for (c in 0:11) {
  print(c)
  x <- df[1:n_top, paste0('gene_NMF', c)]  %>% unlist(use.names = FALSE)
  list.genes[[paste0('NMF', c)]] <- x
}

genes <- c(list.genes$NMF0, list.genes$NMF1, list.genes$NMF2, list.genes$NMF3, list.genes$NMF4,
           list.genes$NMF5, list.genes$NMF6, list.genes$NMF7, list.genes$NMF8, list.genes$NMF9,
           list.genes$NMF10, list.genes$NMF11)
d.plot <- df.corr[genes,]

p <- pheatmap(t(d.plot), cluster_rows = FALSE, cluster_cols = FALSE,
              gaps_row=seq(0,11), gaps_col = seq(0, 12*n_top-1, n_top),
              color =ocean.balance(100), breaks = seq(-1, 1, by=0.02), 
              labels_col=genes, labels_row = row.labels)
p

pdf('../scanpy/220438_NMFprojection/graph/pheatmap.corr.top10.pdf', width = 16, height = 3)
p
dev.off()

sessionInfo()
