library(tidyverse)
library(ggbeeswarm)

setwd('~/autoimmune_10x/')

project.name <- "GSE132338"


df.metadata <- read_csv(paste0('data/', project.name,'/', project.name,'_metadata.csv'))
df.output <- read_csv(paste0('output/', project.name, '/',  project.name, '_queryL2_Reference_Mapping.csv'))

df.output <- df.output %>% 
  left_join(df.metadata, by = "sample") %>%
  group_by(sample) %>%
  mutate(freq = n / sum(n))

p <- ggplot(df.output, aes(x=clusterL2, y=freq, color=disease), scale = "width") + 
  geom_quasirandom(dodge.width=1) +
  theme(axis.text.x = element_text(angle = 90))
p  

ggsave(paste0('output/', project.name,'/', project.name, '_queryL2_ggbeeswarm.pdf'), width = 8, height = 4)


# sample number (per disease)
df.metadata %>% 
  filter(sample %in% df.output$sample) %>%
  group_by(disease) %>%
  summarise(n())

# CD4T number
sum(df.output$n)


# set value 0 in the cluster in the samples was not detected
df.output <- df.output %>% 
  pivot_wider(values_from = freq, names_from = clusterL2,
              id_cols = c(sample, disease), values_fill = 0) %>%
  pivot_longer(!c(sample,disease), names_to = "clusterL2", values_to = "freq")

# change disease name depending on the dataset
dis <- "Sarcoidosis"
df.ttest <- df.output  %>% 
  group_by(clusterL2, disease) %>% 
  summarise(freq = list(freq)) %>% 
  spread(disease, freq) %>% 
  group_by(clusterL2) %>% 
  mutate(p_value = t.test(unlist(.data[[dis]]), unlist(HC))$p.value,
         t_value = t.test(unlist(.data[[dis]]), unlist(HC))$statistic)

df.ttest <- df.ttest %>% select(c(clusterL2, p_value, t_value))
df.ttest['p.adj'] <- p.adjust(df.ttest$p_value)
write_csv(df.ttest, paste0('output/', project.name,'/', project.name, '_queryL2_ttest.csv'))
