---
title: "heatmap_jia"
output: html_document
---

```{read in data}
P0_2 <- read.csv("/project/bf528/project_2/data/samples/P0_2/genes.fpkm_tracking", header=TRUE, sep="\t") 
P0_1 <- read.csv("/projectnb/bf528/users/group_5/project_2/evie/genes.fpkm_tracking", header=TRUE, sep="\t") 
P4_1 <- read.csv("/project/bf528/project_2/data/samples/P4_1/genes.fpkm_tracking", header=TRUE, sep="\t") 
P4_2 <- read.csv("/project/bf528/project_2/data/samples/P4_2/genes.fpkm_tracking", header=TRUE, sep="\t")
P7_1 <- read.csv("/project/bf528/project_2/data/samples/P7_1/genes.fpkm_tracking", header=TRUE, sep="\t")
P7_2 <- read.csv("/project/bf528/project_2/data/samples/P7_2/genes.fpkm_tracking", header=TRUE, sep="\t") 

Ad1_ftable <-read.csv("/project/bf528/project_2/data/samples/Ad_1/genes.fpkm_tracking", header=TRUE, sep="\t")
Ad2_ftable <- read.csv("/project/bf528/project_2/data/samples/Ad_2/genes.fpkm_tracking", header=TRUE, sep="\t")
```


```{get top 1000 genes}
selectg=gene_exp[-which(gene_exp$log2.fold_change. == Inf),]
selectg=selectg[-which(selectg$log2.fold_change. == -Inf),]
selectg=selectg[-which(selectg$significant == 'no'),]
selectg=selectg[c('gene', 'log2.fold_change.','gene_id')]
selectg$abs_diff=abs(selectg$log2.fold_change.)
selectg=selectg[order(selectg$abs_diff, decreasing= T), ]
topgene=selectg[1:1000,]
topg=topgene[c('gene','gene_id')]
```

```{concatenate datasets}
library(tidyverse)
library(tibble)
#FPKM values
P01_fpkm <- P0_1 %>% filter(gene_short_name %in% topg$gene) %>% select(gene_short_name,FPKM) 
P01_temp <- P01_fpkm[!duplicated(P01_fpkm$gene_short_name), ] %>% rename(P0_1 = FPKM)
P02_fpkm <- P0_2 %>% filter(gene_short_name %in% topg$gene) %>% select(gene_short_name,FPKM) 
P02_temp <- P02_fpkm[!duplicated(P02_fpkm$gene_short_name), ] %>% rename(P0_2 = FPKM)
P42_fpkm <- P4_2 %>% filter(gene_short_name %in% topg$gene) %>% select(gene_short_name,FPKM) 
P42_temp <- P42_fpkm[!duplicated(P42_fpkm$gene_short_name), ] %>% rename(P4_2 = FPKM)
P41_fpkm <- P4_1 %>% filter(gene_short_name %in% topg$gene) %>% select(gene_short_name,FPKM) 
P41_temp <- P41_fpkm[!duplicated(P02_fpkm$gene_short_name), ] %>% rename(P4_1 = FPKM)
P71_fpkm <- P7_1 %>% filter(gene_short_name %in% topg$gene) %>% select(gene_short_name,FPKM) 
P71_temp <- P71_fpkm[!duplicated(P71_fpkm$gene_short_name), ] %>% rename(P7_1 = FPKM)
P72_fpkm <- P7_2 %>% filter(gene_short_name %in% topg$gene) %>% select(gene_short_name,FPKM) 
P72_temp <- P72_fpkm[!duplicated(P71_fpkm$gene_short_name), ] %>% rename(P7_2 = FPKM)
ad1_fpkm <- Ad1_ftable %>% filter(gene_short_name %in% topg$gene) %>% select(gene_short_name,FPKM) 
ad1_temp <- ad1_fpkm[!duplicated(ad1_fpkm$gene_short_name), ] %>% rename(Ad_1 = FPKM)
ad2_fpkm <- Ad2_ftable %>% filter(gene_short_name %in% topg$gene) %>% select(gene_short_name,FPKM) 
ad2_temp <- ad2_fpkm[!duplicated(ad2_fpkm$gene_short_name), ] %>% rename(Ad_2 = FPKM)
#concatenated table
concat04 <- inner_join(P01_temp, P02_temp, by="gene_short_name") %>% inner_join(., P41_temp, by="gene_short_name") %>% inner_join(., P42_temp, by="gene_short_name") %>% inner_join(., P71_temp, by="gene_short_name") %>% inner_join(., P72_temp, by="gene_short_name") %>% inner_join(., ad1_temp, by="gene_short_name") %>% inner_join(., ad2_temp, by="gene_short_name")
concat_temp <- concat04[!duplicated(concat04$gene_short_name), ] 
rownames(concat_temp) <- NULL
concat_fin <- concat_temp %>% column_to_rownames("gene_short_name")
```

```{visualization}
sampleseven=as.matrix(concat_fin)
heatmap(sampleseven)
```
