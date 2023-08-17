library(tidyverse)
library(DESeq2)
library(EnhancedVolcano)
library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)

setwd("~/media32TB/bioinformatics/autoimmune_10x/")

# threshold for DESeq2
smallestGroupSize <- 3
th.exp <- 10
# threshold for DEGs
th.padj = 0.1
th.lfc = 1

df.cli <- read_tsv('./data/E-GEAD-397_ImmuNexUT/clinical_diagnosis_age_sex_v2.txt')

vec.disease <- df.cli %>% 
  dplyr::select(disease) %>% 
  unique() %>% 
  filter(! disease %in% c("HC", "BD", "AOSD", "AAV", "SjS")) %>% 
  as.vector()
vec.disease <- vec.disease$disease

vec.cell <- c("Naive_CD4", "Mem_CD4", "Fr_I_nTreg", "Fr_II_eTreg", "Fr_III_T", "Tfh", 'Th1', "Th2", 'Th17')

for (dis in vec.disease) {

list.geneList <- list()
for (cell in vec.cell) {
df <- read_tsv(paste0('./data/E-GEAD-397_ImmuNexUT/',cell ,'_count.txt'))
df <- df %>% 
  dplyr::select(!Gene_id) %>%
  group_by(Gene_name) %>%
  summarise(across(everything(), sum)) %>%
  column_to_rownames(var = "Gene_name")

df.cli.query <- df.cli %>% 
  filter(disease %in% c('HC', dis)) %>%
  filter(id %in% colnames(df))

dds <- DESeqDataSetFromMatrix(countData=df[df.cli.query$id], 
                              colData=df.cli.query, 
                              design=~disease)
dds$condition <- relevel(dds$disease, ref = "HC")

# filter genes

keep <- rowSums(counts(dds) >= th.exp) >= smallestGroupSize
dds <- dds[keep,]

dds <- DESeq(dds)
# res <- results(dds, lfcThreshold = th.lfc, alpha = th.padj)
res <- results(dds)
summary(res)

res.tidy <- results(dds, tidy=TRUE)
write_csv(res.tidy, paste0('./output/Immunexut/', dis, '_', cell, 'DESeq2.csv'))

EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'padj',
                ylab=bquote(~-Log[10] ~ italic(Padj)),
                pCutoff = th.padj,
                FCcutoff = th.lfc,
                title = paste0(dis, ' - ', cell),
                subtitle = paste0('Up ', dim(res.tidy %>% filter((padj < th.padj)&(log2FoldChange>th.lfc)))[1], 'genes ',
                                  'Down ', dim(res.tidy %>% filter((padj < th.padj)&(log2FoldChange < -th.lfc)))[1], 'genes')
)
# ggsave(paste0('./output/Immunexut/', dis, '_', cell, '.volcano.pdf'), width = 8, height = 8)
ggsave(paste0('./output/Immunexut/', dis, '_', cell, '.volcano.png'), width = 8, height = 8)

gL <- res.tidy[c('row', 'log2FoldChange')] %>% 
  # column_to_rownames(var = "row") %>% 
  arrange(desc(log2FoldChange)) %>%
  as.data.frame()
eg <- bitr(gL$row, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
gL <- gL %>% left_join(eg, by=c('row'='SYMBOL')) %>% drop_na() %>% 
  column_to_rownames(var = "ENTREZID")
geneList <- gL$log2FoldChange
names(geneList) <- row.names(gL)

list.geneList[[cell]] <- geneList
}
res <- compareCluster(list.geneList, fun="gsePathway")
dotplot(res, size="NES")
ggsave(paste0('./output/Immunexut/', dis, '.gseReactome.dotplot.pdf'), width = 8, height = 8)
}
