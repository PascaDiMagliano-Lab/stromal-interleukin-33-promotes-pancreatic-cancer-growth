##### DGE ####
# Setup -------------------------------------------------------------------

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
dir.create("outputs", recursive = T, showWarnings = F)

suppressPackageStartupMessages({
  library(DESeq2)
  library(SummarizedExperiment)
  library(limma)
  library(tidyverse)
  library(magrittr)
  source("utils.R")
})


# Loading data ------------------------------------------------------------

counts <- read.table("data/gene_expected_count.annot.txt", sep = "\t", header = T, quote = "") %>%
  `colnames<-`(gsub("^X", "", colnames(.)))
gene_symbols <- counts %>% dplyr::select(gene_id,external_gene_name) %>% column_to_rownames('gene_id')
counts <- counts %>%
  dplyr::select(-c(entrezgene_id, external_gene_name, description)) %>%
  column_to_rownames('gene_id')
metadata <- read.csv("data/samplesheet.csv") %>%
  mutate(sample = gsub("\\-","\\.", sample)) %>%
  column_to_rownames('sample') %>%
  mutate(IL33.Status = as.factor(IL33.Status))
metadata <- metadata[colnames(counts),]
metadata$IL33.Status <- factor(metadata$IL33.Status, levels = c('WT', 'KO'))
metadata$Treatment <- factor(metadata$Treatment, levels = c('NT', 'CM'))

# All samples DESeq object ------------------------------------------------

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = metadata[colnames(counts),],
                              design = ~ 1)
dds <- DESeq(dds)
saveRDS(dds,"outputs/all_samples_dds.rds")

# CM vs. NT (Il33 WT) -----------------------------------------------------

## Data Setup
wt_metadata <- metadata %>%
  filter(IL33.Status == "WT") %>%
  mutate_all(as.factor)
wt_metadata$Treatment <- relevel(wt_metadata$Treatment, ref = "NT")
wt_metadata <- wt_metadata[order(wt_metadata$Treatment),]
wt_counts <- counts[,rownames(wt_metadata)]

## Filtering out genes with low count
wt_counts <- wt_counts[rowSums(wt_counts) != 0,]

## Create DESeq object
wt_design <- model.matrix( ~ Batch + Treatment, data = wt_metadata)
# qr(wt_design)$rank == dim(wt_design)[2]

# creating DESeq object
wt_ddsMat <- DESeqDataSetFromMatrix(countData = wt_counts,
                                    colData = wt_metadata,
                                    design = wt_design)
## optional filtration
wt_ddsMat <- estimateSizeFactors(wt_ddsMat)
idx <- rowSums(counts(wt_ddsMat, normalized=TRUE) >= 10 ) >= 4
wt_ddsMat <- wt_ddsMat[idx,]
wt_dds <- DESeq(wt_ddsMat)
saveRDS(wt_dds, "outputs/treatmentCM_dge_(Il33wt).rds")

lfcShrink(wt_dds, coef = 'TreatmentCM', type = 'ashr') %>%
  data.frame() %>%
  merge(gene_symbols, by = 0) %>%
  filter(!is.na(padj)) %>%
  rename(gene_id = Row.names) %>%
  mutate(status = ifelse(log2FoldChange > 0.5 & padj < 0.05, 1,
                         ifelse(log2FoldChange < -0.5 & padj < 0.05, -1, 0))) %>%
  write.csv("outputs/treatmentCM_dge_(Il33wt).csv", row.names = F, quote = F)

# CM vs. NT (Il33 KO) -----------------------------------------------------

## Data Setup
KO_metadata <- metadata %>%
  filter(IL33.Status == "KO") %>%
  mutate_all(as.factor)
KO_metadata$Treatment <- relevel(KO_metadata$Treatment, ref = "NT")
KO_count <- counts[,rownames(KO_metadata)]

## Filtering out genes with low count
KO_count <- KO_count[rowSums(KO_count) != 0,]

## Create DESeq object
KO_design <- model.matrix( ~ Batch + Treatment, data = KO_metadata)
# qr(KO_design)$rank == dim(KO_design)[2]

# creating DESeq object
KO_ddsMat <- DESeqDataSetFromMatrix(countData = KO_count,
                                    colData = KO_metadata,
                                    design = KO_design)
## Optional Filtration
KO_ddsMat <- estimateSizeFactors(KO_ddsMat)
idx <- rowSums(counts(KO_ddsMat, normalized=TRUE) >= 10 ) >= 4
KO_ddsMat <- KO_ddsMat[idx,]
KO_dds <- DESeq(KO_ddsMat)
saveRDS(KO_dds, "outputs/treatmentCM_dge_(Il33ko).rds")

lfcShrink(KO_dds, coef = 'TreatmentCM', type = 'ashr') %>%
  data.frame() %>%
  merge(gene_symbols, by = 0) %>%
  filter(!is.na(padj)) %>%
  rename(gene_id = Row.names) %>%
  mutate(status = ifelse(log2FoldChange > 0.5 & padj < 0.05, 1,
                         ifelse(log2FoldChange < -0.5 & padj < 0.05, -1, 0))) %>%
  write.csv("outputs/treatmentCM_dge_(Il33ko).csv", row.names = F, quote = F)

# Il33 KO. vs WT (NT) -----------------------------------------------------

## Data Setup
nt_metadata <- metadata %>%
  filter(Treatment == "NT") %>%
  mutate_all(as.factor)
nt_metadata$IL33.Status <- relevel(nt_metadata$IL33.Status, ref = "WT")
nt_count <- counts[,rownames(nt_metadata)]

## Filtering out genes with low count
nt_count <- nt_count[rowSums(nt_count) != 0,]

## Create DESeq object
nt_design <- model.matrix( ~ Batch + IL33.Status, data = nt_metadata)
# qr(nt_design)$rank == dim(nt_design)[2]

# creating DESeq object
nt_ddsMat <- DESeqDataSetFromMatrix(countData = nt_count,
                                    colData = nt_metadata,
                                    design = nt_design)
## Optional Filtration
nt_ddsMat <- estimateSizeFactors(nt_ddsMat)
idx <- rowSums(counts(nt_ddsMat, normalized=TRUE) >= 10 ) >= 4
nt_ddsMat <- nt_ddsMat[idx,]
nt_dds <- DESeq(nt_ddsMat)
saveRDS(nt_dds, "outputs/Il33ko_dge_(NT).rds")

lfcShrink(nt_dds, coef = 'IL33.StatusKO', type = 'ashr') %>%
  data.frame() %>%
  merge(gene_symbols, by = 0) %>%
  filter(!is.na(padj)) %>%
  rename(gene_id = Row.names) %>%
  mutate(status = ifelse(log2FoldChange > 0.5 & padj < 0.05, 1,
                         ifelse(log2FoldChange < -0.5 & padj < 0.05, -1, 0))) %>%
  write.csv("outputs/Il33ko_dge_(NT).csv", row.names = F, quote = F)


# Il33 KO. vs WT (CM) -----------------------------------------------------

cm_metadata <- metadata %>%
  filter(Treatment == "CM") %>%
  mutate_all(as.factor)
cm_metadata$IL33.Status <- relevel(cm_metadata$IL33.Status, ref = "WT")
cm_count <- counts[,rownames(cm_metadata)]

## Filtering out genes with low count
cm_count <- cm_count[rowSums(cm_count) != 0,]

## Create DESeq object
cm_design <- model.matrix( ~ Batch + IL33.Status, data = cm_metadata)
#qr(cm_design)$rank == dim(cm_design)[2]

# creating DESeq object
cm_ddsMat <- DESeqDataSetFromMatrix(countData = cm_count,
                                    colData = cm_metadata,
                                    design = cm_design)
## Optional Filtration
cm_ddsMat <- estimateSizeFactors(cm_ddsMat)
idx <- rowSums(counts(cm_ddsMat, normalized=TRUE) >= 10 ) >= 4
cm_ddsMat <- cm_ddsMat[idx,]
cm_dds <- DESeq(cm_ddsMat)
saveRDS(cm_dds, "outputs/Il33ko_dge_(CM).rds")

lfcShrink(cm_dds, coef = 'IL33.StatusKO', type = 'ashr') %>%
  data.frame() %>%
  merge(gene_symbols, by = 0) %>%
  filter(!is.na(padj)) %>%
  rename(gene_id = Row.names) %>%
  mutate(status = ifelse(log2FoldChange > 0.5 & padj < 0.05, 1,
                         ifelse(log2FoldChange < -0.5 & padj < 0.05, -1, 0))) %>%
  write.csv("outputs/Il33ko_dge_(CM).csv", row.names = F, quote = F)




##### Enrichment Analysis ####
# Setup -------------------------------------------------------------------

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
dir.create("outputs", recursive = T, showWarnings = F)

suppressPackageStartupMessages({
  library(DESeq2)
  library(SummarizedExperiment)
  library(limma)
  library(tidyverse)
  library(magrittr)
  library(fgsea)
  library(msigdbr)
  library(clusterProfiler)
})

csvs <- lapply(list.files("outputs/", pattern = '.csv', full.names = T), read.csv) %>%
  `names<-`(gsub(".csv","",list.files("outputs/", pattern = '.csv')))

# Run GSEA ----------------------------------------------------------------
hallmarks_pathways <- msigdbr(species = "Mus musculus", category = 'H', subcategory = NULL)
hallmarks_pathways <- split(x = hallmarks_pathways$ensembl_gene, f = hallmarks_pathways$gs_name)

lapply(names(csvs), function(x){
  csv <- csvs[[x]]
  rank <- csv %>%
    pull(log2FoldChange) %>%
    `names<-`(csv$gene_id) %>% 
    sort()
  
  gsea <- fgsea(pathways = hallmarks_pathways,
                stats = rank)
  
  saveRDS(gsea, file = paste0("outputs/",gsub("dge","gsea",x),".rds"))
  write.csv(x = gsea[,-8], file = paste0("outputs/",gsub("dge","gsea",x),".csv"), quote = FALSE, row.names = F)
})

#### Plots - Figure S4E ####
## Il33 KO (CM)
il33_ko_CM_dge <- read.csv("outputs/Il33ko_dge_(CM).csv")
il33_ko_CM_rank <- il33_ko_CM_dge %>%
  pull(log2FoldChange) %>%
  `names<-`(il33_ko_CM_dge$gene_id) %>%
  sort()

## Il33 KO (NT)
il33_ko_NT_dge <- read.csv("outputs/Il33ko_dge_(NT).csv")
il33_ko_NT_rank <- il33_ko_NT_dge %>%
  pull(log2FoldChange) %>%
  `names<-`(il33_ko_NT_dge$gene_id) %>%
  sort()

plotEnrichment(hallmarks_pathways[["HALLMARK_TGF_BETA_SIGNALING"]],
               il33_ko_CM_rank) + labs(title="HALLMARK_TGF_BETA_SIGNALING - Il33 KO (CM)")

plotEnrichment(hallmarks_pathways[["HALLMARK_TNFA_SIGNALING_VIA_NFKB"]],
               il33_ko_CM_rank) + labs(title="HALLMARK_TNFA_SIGNALING_VIA_NFKB - Il33 KO (CM)")

plotEnrichment(hallmarks_pathways[["HALLMARK_TGF_BETA_SIGNALING"]],
               il33_ko_NT_rank) + labs(title="HALLMARK_TGF_BETA_SIGNALING - Il33 KO (NT)")

plotEnrichment(hallmarks_pathways[["HALLMARK_TNFA_SIGNALING_VIA_NFKB"]],
               il33_ko_NT_rank) + labs(title="HALLMARK_TNFA_SIGNALING_VIA_NFKB - Il33 KO (NT)")