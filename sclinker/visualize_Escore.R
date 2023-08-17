library(tidyverse)

f_traits <- '../data/ldsc/sumstats.v5.csv'
f_annot <- './gene_scores/categories.txt'
base_her <- 'heritability/Roadmap_U_ABC_h2'
suffix <- '_merged'

setwd('/home/yyasumizu/media32TB/bioinformatics/autoimmune_10x/sclinker/CD4T.NMF.W/')
col.labels <- c('NMF0 Cytotoxic', 'NMF1 Treg', 'NMF2 Th17', 'NMF3 Naiveness', 
                'NMF4 Act', 'NMF5 Th2', 'NMF6 Tfh', 'NMF7 IFN', 'NMF8 Cent. Mem.',
                'NMF9 Thymic Emi.', 'NMF10 Tissue', 'NMF11 Th1')

df.traits <- read.csv(f_traits)

df.traits <- df.traits %>% 
  mutate(Name = str_replace(df.traits$File, '.sumstats', ''))

df.annot <- read.csv(f_annot, header=FALSE, col.names=c("Annotations"))

for (i in 1:nrow(df.traits)) {
  f <- paste(base_her, '/', df.traits[i,'Name'], suffix, '.results', sep='')
  d <- read.table(f, header=TRUE) 
  d <- d %>% mutate(Annotations = str_replace(d$Category, "L2_0", "")) %>%
    filter(Annotations %in% df.annot$Annotations)
  d['Escore'] <- (d$Enrichment - d[d$Annotations=='ALL', 'Enrichment'])
  d <- d %>% mutate(Escore = if_else(Escore < 0, 0, Escore)) %>%
    rename(pEscore = Enrichment_p) %>%
    select(c('Annotations', 'Escore', 'pEscore')) %>%
    filter(Annotations != "ALL")
  d['Trait'] <- df.traits[i,'Name']
  
  if (i == 1){
    df <- d
  } else {
    df <- rbind(df, d)
  }
}

df <- df %>% left_join(df.traits, by=c("Trait" = "Name")) %>% 
  mutate(pEscore = if_else(Escore == 0, 1, pEscore))
df['qEscore'] <- p.adjust(df$pEscore, method = "BH", n = length(df$pEscore))
df['-log10(qEscore)'] <- -log10(df$qEscore)

cols <- paste0('NMF', seq(0,11))

####### Autoimmune
df.plot <- df %>% filter(Autoimmune == 1)
df.plot$Annotations <- factor(df.plot$Annotations, levels = cols) %>%
  fct_recode('NMF0 Cytotoxic' = 'NMF0', 'NMF1 Treg' = 'NMF1', 'NMF2 Th17' = 'NMF2', 
             'NMF3 Naiveness' = 'NMF3', 'NMF4 Act' = 'NMF4', 'NMF5 Th2' = 'NMF5',
             'NMF6 Tfh' = 'NMF6', 'NMF7 IFN' = 'NMF7', 'NMF8 Cent. Mem.' = 'NMF8',
             'NMF9 Thymic Emi.' = 'NMF9', 'NMF10 Tissue' = 'NMF10', 'NMF11 Th1' = 'NMF11')

# mat <- df.plot %>% pivot_wider(names_from = Annotations, values_from = Escore, id_cols = c("Trait")) %>%
#   column_to_rownames(var = "Trait")
# d <- dist(t(scale(t(mat))))
# hc <- hclust(d, "ave")
# # plot(hc)
# rows <- row.names(mat)[hc$order]
order.row <- df.plot %>% group_by(Trait) %>% filter(Escore == max(Escore)) %>% arrange(Annotations) %>% select(Trait)
# order.row <- df.plot %>% group_by(Trait) %>% filter(`-log10(qEscore)`== max(`-log10(qEscore)`)) %>% arrange(Annotations)
df.plot$Trait <- factor(df.plot$Trait, levels = rev(order.row$Trait))
df.plot <- df.plot %>% mutate(Escore = if_else(Escore < 30, Escore, 30))

ggplot(df.plot, aes(Annotations, Trait, size = -log10(qEscore), colour = Escore)) +
  geom_point() +
  scale_color_gradient(low="white", high="red", limit=c(0,30)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave('img/AID_arrange.pdf', width = 8, height = 5)

####### Inflammatory
df.plot <- df %>% filter(Inflammatory == 1)
df.plot$Annotations <- factor(df.plot$Annotations, levels = cols) %>%
  fct_recode('NMF0 Cytotoxic' = 'NMF0', 'NMF1 Treg' = 'NMF1', 'NMF2 Th17' = 'NMF2', 
             'NMF3 Naiveness' = 'NMF3', 'NMF4 Act' = 'NMF4', 'NMF5 Th2' = 'NMF5',
             'NMF6 Tfh' = 'NMF6', 'NMF7 IFN' = 'NMF7', 'NMF8 Cent. Mem.' = 'NMF8',
             'NMF9 Thymic Emi.' = 'NMF9', 'NMF10 Tissue' = 'NMF10', 'NMF11 Th1' = 'NMF11')

# mat <- df.plot %>% pivot_wider(names_from = Annotations, values_from = Escore, id_cols = c("Trait")) %>%
#   column_to_rownames(var = "Trait")
# d <- dist(t(scale(t(mat))))
# hc <- hclust(d, "ave")
# # plot(hc)
# rows <- row.names(mat)[hc$order]
order.row <- df.plot %>% group_by(Trait) %>% filter(Escore == max(Escore)) %>% 
  arrange(Annotations)
df.plot$Trait <- factor(df.plot$Trait, levels = rev(order.row$Trait))


ggplot(df.plot, aes(Annotations, Trait, size = -log10(qEscore), colour = Escore)) +
  geom_point() +
  scale_color_gradient(low="white", high="red") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave('img/Inflamamatory.pdf', width = 8, height = 3)

####### COVID
df.plot <- df %>% filter(COVID == 1)
df.plot$Annotations <- factor(df.plot$Annotations, levels = cols) %>%
  fct_recode('NMF0 Cytotoxic' = 'NMF0', 'NMF1 Treg' = 'NMF1', 'NMF2 Th17' = 'NMF2', 
             'NMF3 Naiveness' = 'NMF3', 'NMF4 Act' = 'NMF4', 'NMF5 Th2' = 'NMF5',
             'NMF6 Tfh' = 'NMF6', 'NMF7 IFN' = 'NMF7', 'NMF8 Cent. Mem.' = 'NMF8',
             'NMF9 Thymic Emi.' = 'NMF9', 'NMF10 Tissue' = 'NMF10', 'NMF11 Th1' = 'NMF11')

# mat <- df.plot %>% pivot_wider(names_from = Annotations, values_from = Escore, id_cols = c("Trait")) %>%
#   column_to_rownames(var = "Trait")
# d <- dist(t(scale(t(mat))))
# hc <- hclust(d, "ave")
# # plot(hc)
# rows <- row.names(mat)[hc$order]
order.row <- df.plot %>% group_by(Trait) %>% filter(Escore == max(Escore)) %>% 
  arrange(Annotations)
df.plot$Trait <- factor(df.plot$Trait, levels = rev(order.row$Trait))

ggplot(df.plot, aes(Annotations, Trait, size = -log10(qEscore), colour = Escore)) +
  geom_point() +
  scale_color_gradient(low="white", high="red") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave('img/COVID.pdf', width = 8, height = 3)


####### Blood
df.plot <- df %>% filter(Blood == 1)
df.plot$Annotations <- factor(df.plot$Annotations, levels = cols) %>%
  fct_recode('NMF0 Cytotoxic' = 'NMF0', 'NMF1 Treg' = 'NMF1', 'NMF2 Th17' = 'NMF2', 
             'NMF3 Naiveness' = 'NMF3', 'NMF4 Act' = 'NMF4', 'NMF5 Th2' = 'NMF5',
             'NMF6 Tfh' = 'NMF6', 'NMF7 IFN' = 'NMF7', 'NMF8 Cent. Mem.' = 'NMF8',
             'NMF9 Thymic Emi.' = 'NMF9', 'NMF10 Tissue' = 'NMF10', 'NMF11 Th1' = 'NMF11')

ggplot(df.plot, aes(Annotations, Trait, size = -log10(qEscore), colour = Escore)) +
  geom_point() +
  scale_color_gradient(low="white", high="red") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave('img/blood.pdf', width = 8, height = 4)

####### Brain
df.plot <- df %>% filter(Brain == 1)
df.plot$Annotations <- factor(df.plot$Annotations, levels = cols)

ggplot(df.plot, aes(Annotations, Trait, size = -log10(qEscore), colour = Escore)) +
  geom_point() +
  scale_color_gradient(low="white", high="red") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


####### Neurodegenerative
df.plot <- df %>% filter(Neurodegenerative == 1)
df.plot$Annotations <- factor(df.plot$Annotations, levels = cols)

ggplot(df.plot, aes(Annotations, Trait, size = -log10(qEscore), colour = Escore)) +
  geom_point() +
  scale_color_gradient(low="white", high="red") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


######## COVID
df.plot <- df %>% filter(grepl('COVID', Trait))
df.plot$Annotations <- factor(df.plot$Annotations, levels = cols)

ggplot(df.plot, aes(Annotations, Trait, size = -log10(qEscore), colour = Escore)) +
  geom_point() +
  scale_color_gradient(low="white", high="red") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

