library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(pals)
library(ggbeeswarm)


setwd("~/media32TB/bioinformatics/autoimmune_10x/NMFprojection/")
row.labels <- c('NMF0 Cytotoxic', 'NMF1 Treg', 'NMF2 Th17', 'NMF3 Naiveness', 
                'NMF4 Act', 'NMF5 Th2', 'NMF6 Tfh', 'NMF7 IFN', 'NMF8 Cent. Mem.',
                'NMF9 Thymic Emi.', 'NMF10 Tissue', 'NMF11 Th1')
# row.labels <- rev(row.labels)

###############
df.proj <- read.csv("DRA008294_iTreg_CD28_projection.csv", row.names = 1)
df.evar <- read.csv("DRA008294_iTreg_CD28_ExplainedVariance.csv", row.names = 1)

# raw value
anno = row_anno_barplot(
  df.evar$ExplainedVariance[1:dim(df.evar)[1]-1],
  border = FALSE, bar_width = 0.9, 
  gp = gpar(fill = cividis(10)[4], lwd = 0),
  width = unit(3, "cm"))
Heatmap(df.proj, name=" ",
        row_order = rownames(df.proj), column_order = colnames(df.proj),
        right_annotation = rowAnnotation(EVar = anno, simple_anno_size = unit(4, "cm")),
        col = cividis(40), row_labels = row.labels, 
        rect_gp = gpar(col = "white", lwd = 2),
        width = ncol(df.proj)*unit(5, "mm"), 
        height = nrow(df.proj)*unit(5, "mm")
)

# scaled value
anno = row_anno_barplot(
  df.evar$ExplainedVariance[1:dim(df.evar)[1]-1],
  border = FALSE, bar_width = 0.9, 
  gp = gpar(fill = ocean.balance(10)[5], lwd = 0),
  width = unit(3, "cm"))
p <- Heatmap(t(scale(t(df.proj))), name=" ",
             row_order = rownames(df.proj), column_order = colnames(df.proj),
             right_annotation = rowAnnotation(EVar = anno, simple_anno_size = unit(4, "cm")),
             col = colorRamp2(seq(-3, 3, length = 20), ocean.balance(20)),
             row_labels = row.labels, 
             rect_gp = gpar(col = "white", lwd = 2),
             width = ncol(df.proj)*unit(5, "mm"), 
             height = nrow(df.proj)*unit(5, "mm")
)

p

pdf("DRA008294_iTreg_CD28.pdf")
p
dev.off()

df.proj.NMF1 <- df.proj %>% as_tibble(rownames = 'NMF') %>%
  pivot_longer(!NMF, names_to = "sample", values_to = "value") %>%
  mutate(group = str_remove(sample, "_0\\d")) %>%
  filter(NMF == 'NMF_1')



###############
df.proj <- read.csv("DRA008294_iTreg_CD28_selected_projection.csv", row.names = 1)
df.evar <- read.csv("DRA008294_iTreg_CD28_selected_ExplainedVariance.csv", row.names = 1)

# raw value
anno = row_anno_barplot(
  df.evar$ExplainedVariance[1:dim(df.evar)[1]-1],
  border = FALSE, bar_width = 0.9, 
  gp = gpar(fill = cividis(10)[4], lwd = 0),
  width = unit(3, "cm"))
Heatmap(df.proj, name=" ",
        row_order = rownames(df.proj), column_order = colnames(df.proj),
        right_annotation = rowAnnotation(EVar = anno, simple_anno_size = unit(4, "cm")),
        col = cividis(40), row_labels = row.labels, 
        rect_gp = gpar(col = "white", lwd = 2),
        width = ncol(df.proj)*unit(5, "mm"), 
        height = nrow(df.proj)*unit(5, "mm")
)

# scaled value
anno = row_anno_barplot(
  df.evar$ExplainedVariance[1:dim(df.evar)[1]-1],
  border = FALSE, bar_width = 0.9, 
  gp = gpar(fill = ocean.balance(10)[5], lwd = 0),
  width = unit(3, "cm"))
p <- Heatmap(t(scale(t(df.proj))), name=" ",
             row_order = rownames(df.proj), column_order = colnames(df.proj),
             right_annotation = rowAnnotation(EVar = anno, simple_anno_size = unit(4, "cm")),
             col = colorRamp2(seq(-3, 3, length = 20), ocean.balance(20)),
             row_labels = row.labels, 
             rect_gp = gpar(col = "white", lwd = 2),
             width = ncol(df.proj)*unit(5, "mm"), 
             height = nrow(df.proj)*unit(5, "mm")
)

p

pdf("DRA008294_iTreg_CD28_selected.pdf")
p
dev.off()

############
df.proj <- read.csv("DICE_projection.csv", row.names = 1)
df.evar <- read.csv("DICE_ExplainedVariance.csv", row.names = 1)

# raw value
anno = row_anno_barplot(
  df.evar$ExplainedVariance[1:dim(df.evar)[1]-1],
  border = FALSE, bar_width = 0.9, 
  gp = gpar(fill = cividis(10)[4], lwd = 0),
  width = unit(3, "cm"))
p <- Heatmap(df.proj, name=" ",
        row_order = rownames(df.proj), column_order = colnames(df.proj),
        right_annotation = rowAnnotation(EVar = anno, simple_anno_size = unit(4, "cm")),
        col = cividis(40), row_labels = row.labels, 
        width = ncol(df.proj)*unit(0.2, "mm"), 
        height = nrow(df.proj)*unit(5, "mm")
)

pdf("DICE.raw.pdf", width = 15)
p
dev.off()

# scaled value
anno = row_anno_barplot(
  df.evar$ExplainedVariance[1:dim(df.evar)[1]-1],
  border = FALSE, bar_width = 0.9, 
  gp = gpar(fill = ocean.balance(10)[5], lwd = 0),
  width = unit(3, "cm"))
p <- Heatmap(t(scale(t(df.proj))), name=" ",
             row_order = rownames(df.proj), column_order = colnames(df.proj),
             right_annotation = rowAnnotation(EVar = anno, simple_anno_size = unit(4, "cm")),
             # col=ocean.balance(20), 
             col = colorRamp2(seq(-3, 3, length = 20), ocean.balance(20)),
             row_labels = row.labels, 
             # rect_gp = gpar(col = "white", lwd = 2),
             width = ncol(df.proj)*unit(0.2, "mm"), 
             height = nrow(df.proj)*unit(5, "mm"),
             column_split = colnames(df.proj) %>% str_remove("_\\d*")
)

p

pdf("DICE.pdf", width = 20)
p
dev.off()

############
df.proj <- read.csv("STR1.5_Fr1.2.3.5.6_projection.csv", row.names = 1)
df.evar <- read.csv("STR1.5_Fr1.2.3.5.6_ExplainedVariance.csv", row.names = 1)
# col.order <- c("act_Tconv_1", "act_Tconv_2", "iTreg_1", "iTreg_2", "nTreg_1", "nTreg_2", 
#                "Ultra_1", "Ultra_2", "Ultraplus_1", "Ultraplus_2")

# raw value
anno = row_anno_barplot(
  df.evar$ExplainedVariance[1:dim(df.evar)[1]-1],
  border = FALSE, bar_width = 0.9, 
  gp = gpar(fill = cividis(10)[4], lwd = 0),
  width = unit(3, "cm"))
Heatmap(df.proj, name=" ",
        row_order = rownames(df.proj), column_order = colnames(df.proj),
        right_annotation = rowAnnotation(EVar = anno, simple_anno_size = unit(4, "cm")),
        col = cividis(40), row_labels = row.labels, 
        rect_gp = gpar(col = "white", lwd = 2),
        width = ncol(df.proj)*unit(5, "mm"), 
        height = nrow(df.proj)*unit(5, "mm")
)

# scaled value
anno = row_anno_barplot(
  df.evar$ExplainedVariance[1:dim(df.evar)[1]-1],
  border = FALSE, bar_width = 0.9, 
  gp = gpar(fill = ocean.balance(10)[5], lwd = 0),
  width = unit(3, "cm"))
p <- Heatmap(t(scale(t(df.proj))), name=" ",
             row_order = rownames(df.proj), column_order = colnames(df.proj),
             right_annotation = rowAnnotation(EVar = anno, simple_anno_size = unit(4, "cm")),
             col = colorRamp2(seq(-3, 3, length = 20), ocean.balance(20)),
             row_labels = row.labels, 
             rect_gp = gpar(col = "white", lwd = 2),
             width = ncol(df.proj)*unit(5, "mm"), 
             height = nrow(df.proj)*unit(5, "mm"),
             column_split = colnames(df.proj) %>% str_extract("Fr.")
)

p

pdf("STR1.5_Fr1.2.3.5.6.pdf")
p
dev.off()

############# Tph

# df.proj <- read.csv("Tph_projection.csv", row.names = 1)
df.proj <- read_csv("Tph_projection.csv") %>% tibble::column_to_rownames(var = '...1')
df.evar <- read.csv("Tph_ExplainedVariance.csv", row.names = 1)
# col.order <- c("act_Tconv_1", "act_Tconv_2", "iTreg_1", "iTreg_2", "nTreg_1", "nTreg_2", 
#                "Ultra_1", "Ultra_2", "Ultraplus_1", "Ultraplus_2")

# raw value
anno = row_anno_barplot(
  df.evar$ExplainedVariance[1:dim(df.evar)[1]-1],
  border = FALSE, bar_width = 0.9, 
  gp = gpar(fill = cividis(10)[4], lwd = 0),
  width = unit(3, "cm"))
Heatmap(df.proj, name=" ",
        row_order = rownames(df.proj), column_order = colnames(df.proj),
        right_annotation = rowAnnotation(EVar = anno, simple_anno_size = unit(4, "cm")),
        col = cividis(40), row_labels = row.labels, 
        rect_gp = gpar(col = "white", lwd = 2),
        width = ncol(df.proj)*unit(5, "mm"), 
        height = nrow(df.proj)*unit(5, "mm")
)

# scaled value
anno = row_anno_barplot(
  df.evar$ExplainedVariance[1:dim(df.evar)[1]-1],
  border = FALSE, bar_width = 0.9, 
  gp = gpar(fill = ocean.balance(10)[5], lwd = 0),
  width = unit(3, "cm"))
p <- Heatmap(t(scale(t(df.proj))), name=" ",
             row_order = rownames(df.proj), column_order = sort(colnames(df.proj)),
             right_annotation = rowAnnotation(EVar = anno, simple_anno_size = unit(4, "cm")),
             col = colorRamp2(seq(-3, 3, length = 20), ocean.balance(20)),
             row_labels = row.labels, 
             rect_gp = gpar(col = "white", lwd = 2),
             width = ncol(df.proj)*unit(5, "mm"), 
             height = nrow(df.proj)*unit(5, "mm"),
             column_split = colnames(df.proj) %>% str_extract("Fr.")
)

p

pdf("Tph.pdf", width = 15)
p
dev.off()

###
df.proj <- read.csv("MS05-07_projection.csv", row.names = 1)
df.evar <- read.csv("MS05-07_ExplainedVariance.csv", row.names = 1)

# raw value
anno = row_anno_barplot(
  df.evar$ExplainedVariance[1:dim(df.evar)[1]-1],
  border = FALSE, bar_width = 0.9, 
  gp = gpar(fill = cividis(10)[4], lwd = 0),
  width = unit(3, "cm"))
p <- Heatmap(df.proj, name=" ",
        row_order = rownames(df.proj), column_order = colnames(df.proj),
        right_annotation = rowAnnotation(EVar = anno, simple_anno_size = unit(4, "cm")),
        col = cividis(40), row_labels = row.labels, 
        rect_gp = gpar(col = "white", lwd = 2),
        width = ncol(df.proj)*unit(5, "mm"), 
        height = nrow(df.proj)*unit(5, "mm")
)

pdf("MS05-07.raw.pdf")
p
dev.off()

# scaled value
anno = row_anno_barplot(
  df.evar$ExplainedVariance[1:dim(df.evar)[1]-1],
  border = FALSE, bar_width = 0.9, 
  gp = gpar(fill = ocean.balance(10)[5], lwd = 0),
  width = unit(3, "cm"))
p <- Heatmap(t(scale(t(df.proj))), name=" ",
             row_order = rownames(df.proj), column_order = colnames(df.proj),
             right_annotation = rowAnnotation(EVar = anno, simple_anno_size = unit(4, "cm")),
             col = colorRamp2(seq(-3, 3, length = 20), ocean.balance(20)),
             row_labels = row.labels, 
             rect_gp = gpar(col = "white", lwd = 2),
             width = ncol(df.proj)*unit(5, "mm"), 
             height = nrow(df.proj)*unit(5, "mm")
)

p

pdf("MS05-07.pdf")
p
dev.off()

############# NCOMM

df.proj <- read_csv("NCOMMS-19-7936188_bulk_RNAseq_projection.csv") %>% tibble::column_to_rownames(var = '...1')
df.evar <- read.csv("NCOMMS-19-7936188_bulk_RNAseq_ExplainedVariance.csv", row.names = 1)
# col.order <- c("act_Tconv_1", "act_Tconv_2", "iTreg_1", "iTreg_2", "nTreg_1", "nTreg_2", 
#                "Ultra_1", "Ultra_2", "Ultraplus_1", "Ultraplus_2")

# raw value
anno = row_anno_barplot(
  df.evar$ExplainedVariance[1:dim(df.evar)[1]-1],
  border = FALSE, bar_width = 0.9, 
  gp = gpar(fill = cividis(10)[4], lwd = 0),
  width = unit(3, "cm"))
Heatmap(df.proj, name=" ",
        row_order = rownames(df.proj), column_order = colnames(df.proj),
        right_annotation = rowAnnotation(EVar = anno, simple_anno_size = unit(4, "cm")),
        col = cividis(40), row_labels = row.labels, 
        rect_gp = gpar(col = "white", lwd = 2),
        width = ncol(df.proj)*unit(5, "mm"), 
        height = nrow(df.proj)*unit(5, "mm")
)

# scaled value
anno = row_anno_barplot(
  df.evar$ExplainedVariance[1:dim(df.evar)[1]-1],
  border = FALSE, bar_width = 0.9, 
  gp = gpar(fill = ocean.balance(10)[5], lwd = 0),
  width = unit(3, "cm"))
p <- Heatmap(t(scale(t(df.proj))), name=" ",
             row_order = rownames(df.proj), column_order = sort(colnames(df.proj)),
             right_annotation = rowAnnotation(EVar = anno, simple_anno_size = unit(4, "cm")),
             col = colorRamp2(seq(-3, 3, length = 20), ocean.balance(20)),
             row_labels = row.labels, 
             rect_gp = gpar(col = "white", lwd = 2),
             width = ncol(df.proj)*unit(5, "mm"), 
             height = nrow(df.proj)*unit(5, "mm"),
             column_split = colnames(df.proj) %>% str_extract("Fr.")
)

p

pdf("NCOMMS-19-7936188_bulk_RNAseq.pdf", width = 15)
p
dev.off()


############# pTreg

df.proj <- read_csv("PRJNA735530_pTreg_bulkRNAseq_projection.csv") %>% tibble::column_to_rownames(var = '...1')
df.evar <- read.csv("PRJNA735530_pTreg_bulkRNAseq_ExplainedVariance.csv", row.names = 1)
# col.order <- c("act_Tconv_1", "act_Tconv_2", "iTreg_1", "iTreg_2", "nTreg_1", "nTreg_2", 
#                "Ultra_1", "Ultra_2", "Ultraplus_1", "Ultraplus_2")

# raw value
anno = row_anno_barplot(
  df.evar$ExplainedVariance[1:dim(df.evar)[1]-1],
  border = FALSE, bar_width = 0.9, 
  gp = gpar(fill = cividis(10)[4], lwd = 0),
  width = unit(3, "cm"))
Heatmap(df.proj, name=" ",
        row_order = rownames(df.proj), column_order = colnames(df.proj),
        right_annotation = rowAnnotation(EVar = anno, simple_anno_size = unit(4, "cm")),
        col = cividis(40), row_labels = row.labels, 
        rect_gp = gpar(col = "white", lwd = 2),
        width = ncol(df.proj)*unit(5, "mm"), 
        height = nrow(df.proj)*unit(5, "mm")
)

# scaled value
anno = row_anno_barplot(
  df.evar$ExplainedVariance[1:dim(df.evar)[1]-1],
  border = FALSE, bar_width = 0.9, 
  gp = gpar(fill = ocean.balance(10)[5], lwd = 0),
  width = unit(3, "cm"))
p <- Heatmap(t(scale(t(df.proj))), name=" ",
             row_order = rownames(df.proj), column_order = sort(colnames(df.proj)),
             right_annotation = rowAnnotation(EVar = anno, simple_anno_size = unit(4, "cm")),
             col = colorRamp2(seq(-3, 3, length = 20), ocean.balance(20)),
             row_labels = row.labels, 
             rect_gp = gpar(col = "white", lwd = 2),
             width = ncol(df.proj)*unit(5, "mm"), 
             height = nrow(df.proj)*unit(5, "mm"),
             column_split = colnames(df.proj) %>% str_extract("Fr.")
)

p

pdf("PRJNA735530_pTreg.pdf", width = 15)
p
dev.off()



######### immunexut
setwd("~/media32TB/bioinfordf.projics/autoimmune_10x/NMFprojection/")

df.proj <- read.csv("ImmuNexUT/CD4T_CD4T_projection.csv", row.names = 1)
df.evar <- read.csv("ImmuNexUT/CD4T_CD4T_ExplainedVariance.csv", row.names = 1)

# raw value
anno = row_anno_barplot(
  df.evar$ExplainedVariance[1:dim(df.evar)[1]-1],
  border = FALSE, bar_width = 0.9, 
  gp = gpar(fill = cividis(10)[4], lwd = 0),
  width = unit(3, "cm"))
p <- Heatmap(df.proj, name=" ",
             row_order = rownames(df.proj), column_order = colnames(df.proj),
             right_annotation = rowAnnotation(EVar = anno, simple_anno_size = unit(4, "cm")),
             col = cividis(40), row_labels = row.labels, 
             width = ncol(df.proj)*unit(0.05, "mm"), 
             height = nrow(df.proj)*unit(5, "mm")
)

p
pdf("ImmuNexUT/img/CD4T_CD4T.raw.pdf", width = 15)
p
dev.off()

# scaled value
anno = row_anno_barplot(
  df.evar$ExplainedVariance[1:dim(df.evar)[1]-1],
  border = FALSE, bar_width = 0.9, 
  gp = gpar(fill = ocean.balance(10)[5], lwd = 0),
  width = unit(3, "cm"))
p <- Heatmap(t(scale(t(df.proj))), name=" ",
             row_order = rownames(df.proj), column_order = colnames(df.proj),
             right_annotation = rowAnnotation(EVar = anno, simple_anno_size = unit(4, "cm")),
             # col=ocean.balance(20), 
             col = colorRamp2(seq(-3, 3, length = 20), ocean.balance(20)),
             row_labels = row.labels, 
             # rect_gp = gpar(col = "white", lwd = 2),
             width = ncol(df.proj)*unit(0.05, "mm"), 
             height = nrow(df.proj)*unit(5, "mm"),
             column_split = colnames(df.proj) %>% str_remove("NW\\d*_")
)

p

pdf("ImmuNexUT/img/CD4T_CD4T.pdf", width = 15)
p
dev.off()

cols <- paste0('NMF_', seq(0,11))
cell <- "Fr_I_nTreg"
cell <- 'Naive_CD4'
df.stats <- read_csv('ImmuNexUT/stats_cd4T_cd4T.csv') %>% 
  filter(index != "const") %>% 
  filter(cell == !!cell) %>%
  mutate(`-log10q`  = if_else(-log10(qvalue) > 10, 10, -log10(qvalue)))
  # mutate(`-log10q` = -log10(qvalue))
df.stats$NMF <- factor(df.stats$NMF, levels = cols)
levels(df.stats$NMF) <- row.labels
# ggplot(df.stats, aes(index, NMF, colour = Coef., size=`-log10q`)) +
#   geom_point() +
#   scale_color_gradientn(colours = rev(brewer.rdbu(20)), limits = c(-0.6, 0.6)) +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
#   ggtitle(cell) 
# ggsave(paste0('ImmuNexUT/img/CD4T_CD4T_', cell, '.pdf'), width = 5, height = 5)

for (cell in c('Fr_III_T', 'Fr_II_eTreg', 'Fr_I_nTreg', 'Mem_CD4', 'Naive_CD4',
                   'Tfh', 'Th17', 'Th1', 'Th2')) {
  df.stats <- read_csv('ImmuNexUT/stats_cd4T_cd4T.csv') %>% 
    filter(index != "const") %>% 
    filter(cell == !!cell) %>%
    mutate(`-log10q`  = if_else(-log10(qvalue) > 10, 10, -log10(qvalue)))
  # mutate(`-log10q` = -log10(qvalue))
  df.stats$NMF <- factor(df.stats$NMF, levels = cols)
  levels(df.stats$NMF) <- row.labels
  ggplot(df.stats, aes(index, NMF, colour = Coef., size=`-log10q`)) +
    geom_point() +
    scale_color_gradientn(colours = rev(brewer.rdbu(20)), limits = c(-0.6, 0.6)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    ggtitle(cell) 
  ggsave(paste0('ImmuNexUT/img/CD4T_CD4T_', cell, '.pdf'), width = 5, height = 5)
} 

