library(dplyr)
library(tidyverse)
library(ggalluvial)

gene_annotations_file = "gene_annot_single_kallisto.tsv"
RNAseq_file = "RNAseq_quantile_normalized_dat.tsv"

LADs = c("T1-LAD", "T2-LAD", "nonLAD")

ct_replace = c(
  'CardiacMyocytes'='Cardiac\nmyocytes',
  'EarlySomite'='Early somite',
  'H9ESC'='ESCs',
  'ParaxMesoderm'='Paraxial\nmesoderm',
  'DefEctoderm'='Definitive\nectoderm',
  'MidHindgut'='Mid-hindgut',
  'D5Midbrain'='Midbrain',
  'Liver'='Liver progenitors',
  'BorderEctoderm'='Border ectoderm',
  'EndoProgenitor'='Endothelial\nprogenitors',
  'D4Artery'='Artery progenitors',
  'APS'='Anterior primitive\nstreak'
)

cell_types = c("H9ESC", "MidHindgut", "Liver", "APS", "ParaxMesoderm", "EarlySomite", 
               "D4Artery", "EndoProgenitor", "CardiacMyocytes", "Epicardium", 
               "DefEctoderm", "BorderEctoderm", "D5Midbrain")
cell_types_RNAseq = c("H9ESC", "MidHindgut", "ParaxMesoderm", "EarlySomite", 
                      "EndoProgenitor", "CardiacMyocytes", "DefEctoderm")

# LAD bin calls
lb1_dat <- read.table("LADs_per_bin_may2021.tsv", sep = "\t", header = T) %>%
  select(chrom, start, stop, category, cell_type = ct) %>%
  mutate(chrom = gsub("chr", "", chrom)) %>%
  mutate(block = stop/20000)



# Annotate genes in all cell types ####

# Read in gene annotations from gencode v37 annotations (only genes that are listed once)
gene_annot <- read.table(gene_annotations_file, header = T, sep = "\t")

gene_annot_sub <- gene_annot %>%
  mutate(chrom = gsub("chr", "", chr)) %>%
  mutate(block_start = if_else(strand == "+", 
                               ceiling(as.numeric(start)/20000), 
                               ceiling(as.numeric(end)/20000))) %>% # start is TSS
  mutate(block_end = if_else(strand == "+", 
                             ceiling(as.numeric(end)/20000), 
                             ceiling(as.numeric(start)/20000))) %>%
  select(gene = name, chrom, block_start, block_end)

gene_annot_LAD <- gene_annot_sub[rep(seq_len(nrow(gene_annot_sub)), 
                                     each = length(cell_types)), ] %>%
  mutate(cell_type = rep(cell_types, times = nrow(gene_annot_sub)))

gene_annot_LAD <- merge(gene_annot_LAD, 
                        select(lb1_dat,
                               cell_type,
                               chrom, 
                               block_start = block, 
                               category_start = category),
                        by = c("chrom", "block_start", "cell_type"),
                        all.x = T, all.y = F) %>%
  select(-block_start, -block_end)



# Get tpm fold change in ESC vs differentiating cells ####


# Get quantile normalized gene expression in tpm
RNAseq_data <- read.table(RNAseq_file, header = T, sep = "\t")

names(RNAseq_data) <- c("gene_ID", "tpm_CardiacMyocytes", "tpm_H9ESC", "tpm_EarlySomite", "tpm_ParaxMesoderm", "tpm_MidHindgut", "tpm_EndoProgenitor", "tpm_DefEctoderm")

cell_types_RNAseq <- c("CardiacMyocytes", "EndoProgenitor", "EarlySomite","ParaxMesoderm",
                       "MidHindgut", "DefEctoderm", "H9ESC")

RNAseq_data_sub <- RNAseq_data %>%
  # remove genes that are not expressed in any of the cell types
  filter(!(tpm_CardiacMyocytes == 0 &
             tpm_H9ESC == 0 &
             tpm_EarlySomite == 0 &
             tpm_ParaxMesoderm == 0 &
             tpm_MidHindgut == 0 &
             tpm_EndoProgenitor == 0 &
             tpm_DefEctoderm == 0))

# Change tpm of 0 to 0.0001 so you can divide 
RNAseq_data_sub[RNAseq_data_sub == 0] <- 0.0001

for (cell in setdiff(cell_types_RNAseq, "H9ESC")){
  RNAseq_data_sub[,cell] <- 
    RNAseq_data_sub[,names(RNAseq_data_sub)[grepl(cell, names(RNAseq_data_sub))]]/
    RNAseq_data_sub[,"tpm_H9ESC"]
}

# Annotate gene expression change categories
fold_change_cont <- pivot_longer(RNAseq_data_sub,
                                 cols = setdiff(cell_types_RNAseq, "H9ESC"),
                                 names_to = "cell_type",
                                 values_to = "fold_change") %>%
  mutate(cell = gsub("_foldchange", "", cell)) %>%
  select(gene = gene_ID, cell_type, fold_change)
  


# get gene LAD calls ####

gene_annot_sub1 <- gene_annot_sub %>%
  filter(block_start != block_end)

gene_blocks <- tibble(Gene = "",
                      chrom = "",
                      block = numeric())

# get LAD assignments for each block of genes that fall on multiple blocks
for (gene in unique(gene_annot_sub1$gene)){
  
  chrom <- gene_annot_sub1[gene_annot_sub1$gene == gene,"chrom"]
  start <- gene_annot_sub1[gene_annot_sub1$gene == gene,"block_start"]
  end <- gene_annot_sub1[gene_annot_sub1$gene == gene,"block_end"]
  
  gene_block <- tibble(Gene = rep(gene, times = length(start:end)),
                       chrom = chrom,
                       block = start:end)
  
  gene_blocks <- rbind(gene_blocks, gene_block)
  
}

gene_blocks2 <- gene_blocks[rep(seq_len(nrow(gene_blocks)), 
                                each = length(cell_types)), ] %>%
  mutate(cell_type = rep(cell_types, times = nrow(gene_blocks)))

gene_blocks_LAD <- merge(
  gene_blocks2, 
  select(lb1_dat, chrom, block, category, cell_type),
  by = c("chrom", "block", "cell_type"),
  all.x = T, all.y = F) %>%
  filter(!is.na(category))



# 5A: Expression vs LAD change ####

fold_change_cont_LAD <- merge(
  
  # Add ESC LAD assignments
  merge(fold_change_cont,
        select(filter(gene_annot_LAD, cell_type == "H9ESC"), gene, category_start),
        by = "gene",
        all.x = T) %>%
    dplyr::rename(ESC_category = category_start),
  
  # Add all other cell type LAD assignments
  select(filter(gene_annot_LAD, cell_type != "H9ESC"), gene, cell_type, category_start),
  
  by = c("gene", "cell_type"), all.x = T) %>%
  dplyr::rename(cell_category = category_start)

# Narrow down based on change in gene expression and LAD assignment of interest

fold_change_cont_LAD <- fold_change_cont_LAD %>%
  filter(!is.na(ESC_category)) %>% # Remove genes that are not mapped to LADs in ESCs ~3.7%
  mutate(LAD_change = paste(ESC_category, cell_category, sep = "_to_")) %>%
  filter(ESC_category != cell_category) %>%
  mutate(LAD_change_2 = if_else(LAD_change == "nonLAD_to_T1-LAD" |
                                  LAD_change == "nonLAD_to_T2-LAD" |
                                  LAD_change == "T2-LAD_to_T1-LAD",
                                "deactivate", "activate"))

fold_change_cont_LAD$LAD_change <- factor(fold_change_cont_LAD$LAD_change, 
                                          levels = c("T1-LAD_to_nonLAD","T1-LAD_to_T2-LAD", "T2-LAD_to_nonLAD",
                                                     "nonLAD_to_T2-LAD", "T2-LAD_to_T1-LAD", "nonLAD_to_T1-LAD"))


# for genes with only one LAD assignment

fold_change_cont_LAD_1 <- fold_change_cont_LAD

# Keep only genes with one assignment
for (cell in unique(fold_change_cont_LAD_1$cell_type)){
  
  # get genes with more than one assignment for each cell type to remove
  to_remove <- (unique(gene_blocks_LAD %>%
                         filter(cell_type == cell) %>%
                         select(Gene, category)) %>%
                  group_by(Gene) %>% 
                  filter(n() > 1))$Gene
  
  fold_change_cont_LAD_1 <- fold_change_cont_LAD_1 %>%
    filter(!(cell_type == cell & gene %in% to_remove))
  
  
}


graph <- fold_change_cont_LAD_1 %>%
  group_by(cell_type, LAD_change) %>%
  summarise(mean = mean(log(fold_change)), 
            sd = sd(log(fold_change)),
            count = n()) %>%
  mutate(low_sd = mean - sd, high_sd = mean + sd) %>%
  filter(cell_type != "EarlySomite" &
           cell_type != "ParaxMesoderm") %>%
  mutate(cell_type = if_else(cell_type == "CardiacMyocytes", "Cardiac myocytes", cell_type)) %>%
  mutate(cell_type = if_else(cell_type == "DefEctoderm", "Definitive ectoderm", cell_type)) %>%
  mutate(cell_type = if_else(cell_type == "EndoProgenitor", "Endothelial progenitors", cell_type)) %>%
  mutate(cell_type = if_else(cell_type == "MidHindgut", "Mid-hindgut", cell_type))

ggplot(graph, aes(x = LAD_change, y = mean, group = 1)) +
  geom_ribbon(aes(ymin=low_sd, ymax=high_sd), alpha=0.2) +
  facet_wrap(~ cell_type, scales = "free", nrow = 1) +
  geom_point(aes(color = LAD_change), size = 3) +
  scale_color_brewer(palette="Dark2") +
  labs(y = "Log(fold change tpm)") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        strip.text = element_text(size=10),
        axis.ticks = element_line())


# Anova test for each cell type
for (cell in unique(fold_change_cont_LAD_1$cell_type)){
  aov <- aov(log(fold_change) ~ LAD_change, 
             data = filter(fold_change_cont_LAD_1, cell_type == cell))
  pval <- summary(aov)[[1]][1,5]
  print(cell)
  print(pval)
  print(filter(as.data.frame(TukeyHSD(aov)$LAD_change), `p adj` < 0.05))
}



# 5B: Characteristic genes ####

temp <- gene_annot_LAD %>%
  mutate(annotation = if_else((cell_type == "CardiacMyocytes" & 
                                 gene %in% c("TNNT2", "ACTC1", "MYH6", "MYL7", "MYBPC3")) |
                                (cell_type == "DefEctoderm" & 
                                   gene %in% c("OTX2", "BMP4", "PAX6", "SOX1", "FOXJ3")) |
                                (cell_type == "MidHindgut" & 
                                   gene %in% c("CDX2", "FOXA2", "HOXA9", "HOXC5", "SOX17")) |
                                (cell_type == "BorderEctoderm" & 
                                   gene %in% c("PAX3", "DLX5", "ZIC1", "MSX1", "PAX7")) |
                                (cell_type == "D5Midbrain" & 
                                   gene %in% c("PAX2", "PAX5", "EN1", "EN2", "TCF4")) |
                                (cell_type == "D4Artery" & 
                                   gene %in% c("PECAM1", "CD34", "CDH5", "DLL4", "SOX17")) |
                                (cell_type == "Liver" & 
                                   gene %in% c("HNF4A", "TBX3", "APOH", "F9", "FGA")) |
                                ((cell_type == "ParaxMesoderm" | cell_type == "EarlySomite") & 
                                   gene %in% c("CDX2", "TBX6", "MSGN1", "MEOX1", "MSX1")) |
                                (cell_type == "EndoProgenitor" & 
                                   gene %in% c("CD34", "CXCR4", "TEK", "KDR", "PECAM1")) |
                                (cell_type == "H9ESC" & 
                                   gene %in% c("NANOG", "SOX2", "KLF4", "MYC", "ZFP42")),
                              "X", "")
  ) %>%
  filter(cell_type %in% cell_types_RNAseq) %>%
  arrange(factor(cell_type, levels = cell_types_RNAseq))

# list of characteristic genes ordered by cell type
genes_ofinterest <- unique(filter(temp, annotation == "X")$gene)

graph <- temp %>% 
  filter(gene %in% genes_ofinterest) %>%
  dplyr::rename(category = category_start)

# Incorporate gene expression data

RNAseq_data_ofinterest <- pivot_longer(RNAseq_data, 
                                       cols = names(RNAseq_data)[grepl("tpm", names(RNAseq_data))],
                                       names_to = "cell_type",
                                       values_to = "tpm") %>%
  mutate(cell_type = gsub("tpm_", "", cell_type)) %>%
  dplyr::rename(gene = gene_ID) %>%
  filter(gene %in% genes_ofinterest) %>%
  mutate(tpm = if_else(tpm == 0, 0.0001, tpm))

graph1 <- merge(graph, RNAseq_data_ofinterest, 
                by = c("cell_type", "gene"),
                all.x = T) %>%
  select(gene, cell_type, category, annotation, tpm) %>%
  arrange(factor(cell_type, levels = cell_types_RNAseq))

graph1$cell_type <- str_replace_all(graph1$cell_type, ct_replace)
graph1$cell_type <- factor(graph1$cell_type, levels = rev(unique(graph1$cell_type)))

# Order category and genes
graph1$category <- factor(graph1$category, levels = LADs)
graph1$gene <- factor(graph1$gene, levels = genes_ofinterest)

ggplot(graph1, aes(y = cell_type, x = gene, shape = category)) +
  geom_point(aes(color = log10(tpm)), size = 4) +
  scale_color_gradient(low = "lightgrey", high = "goldenrod3", na.value="snow2") +
  geom_text(aes(label = annotation), color = "black") +  
  labs(shape = "Category", color = "Log10(tpm)") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_text(size=12, face = "bold"),
        legend.text = element_text(size=10))


# 5C: Genes across trajectory - Alluvial plot ####

cell_1 = "H9ESC"
cell_2 = "ParaxMesoderm"
cell_3 = "CardiacMyocytes"

lineage = c(cell_1, cell_2, cell_3)

# Get gene LAD calls for each of the 3 cell types
lineage_LADs <- pivot_wider(
  gene_annot_LAD %>%
    select(gene, cell_type, category = category_start) %>%
    filter(!is.na(category)) %>%
    filter(cell_type %in% lineage) %>%
    mutate(cell_order = ifelse(cell_type == cell_1, "c1", "")) %>%
    mutate(cell_order = ifelse(cell_type == cell_2, "c2", cell_order)) %>%
    mutate(cell_order = ifelse(cell_type == cell_3, "c3", cell_order)) %>%
    select(-`cell_type`),
  names_from = "cell_order",
  values_from = "category"
) %>%
  mutate(switch = paste(c1, c2, c3, sep = "_"))

graph <- unique(merge(lineage_LADs,
                      dplyr::rename(
                        as.data.frame(table(lineage_LADs$switch)),
                        switch = Var1),
                      by = "switch") %>%
                  select(-gene, -switch)) %>%
  mutate(Freq = Freq/1000)

# Order by LAD call
for (i in 1:3){
  graph[,i] <- factor(graph[,i], levels = LADs)
}


lineage <- str_replace_all(lineage, ct_replace)

ggplot(graph,
       aes(y = Freq,
           axis1 = c1, axis2 = c2, axis3 = c3)) +
  geom_alluvium(aes(fill = c1),
                width = 0, knot.pos = 0, reverse = FALSE) +
  scale_fill_manual(values = alpha(c("#4B0082", "#9370DB", "lightgrey"), 1)) + # indigo and mediumpurple
  guides(fill = FALSE) +
  geom_stratum(width = 1/8, reverse = FALSE) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), 
            reverse = FALSE, angle = 90, size = 3) +
  scale_x_continuous(breaks = 1:3, labels = lineage) +
  labs(y = "Number of genes (x10k)") +
  theme_classic() +
  theme(axis.text.x = element_text(size=10, color = "black"),
        axis.text.y = element_text(size=10, color = "black"),
        axis.title.y = element_text(size=12, color = "black"),
        axis.title.x = element_blank())


# 5D/E: Heatmaps for tpm of select genes ####

cell_1 = "H9ESC"
cell_2 = "ParaxMesoderm"
cell_3 = "CardiacMyocytes"

lineage = c(cell_1, cell_2, cell_3)

graph <- pivot_longer(RNAseq_data, 
             cols = names(RNAseq_data)[grepl("tpm", names(RNAseq_data))],
             names_to = "cell_type",
             values_to = "tpm") %>%
  mutate(cell_type = gsub("tpm_", "", cell_type)) %>%
  dplyr::rename(gene = gene_ID) %>%
  mutate(tpm = if_else(tpm == 0, 0.0001, tpm)) %>%
  filter(cell_type %in% lineage) %>%
  filter(gene == "SCCPDH") %>% # SCCPDH or TBX20
  mutate(tpm = round(tpm, 1))

graph$cell_type <- factor(graph$cell_type, levels = rev(lineage))

ggplot(graph, aes(x = gene, y = cell_type, color = tpm)) +
  geom_point(size = 15, shape = "square") +
  scale_color_gradient(low = "lightgrey", high = "goldenrod", na.value="snow2") +
  geom_text(aes(label = tpm), color = "black", size = 4) +
  theme_classic() +  
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_blank(),
        legend.title = element_text(size=12, face = "bold"),
        legend.text = element_text(size=10))







