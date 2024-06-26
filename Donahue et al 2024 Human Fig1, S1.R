#This script recreates the human scRNAseq analysis from Figures 1B&C and S1A of Donahue et al. 2024.
#This analysis begins with the Seurat object from the human PDA and adjacent normal datatset originally published in:

#Steele, N.G., Carpenter, E.S., Kemp, S.B. et al. 
#Multimodal mapping of the tumor and peripheral blood immune landscape in human pancreatic cancer. 
#Nat Cancer 1, 1097–1112 (2020). https://doi.org/10.1038/s43018-020-00121-4. 

#Dataset is available on NCBI dbGaP #PHS002071.v1.p1

#Data was processed in line with the Seurat workflow:
#Website: https://satijalab.org/seurat/index.html

#Reference: Hao et al., Integrated analysis of multimodal single-cell data. 
#Cell. 2021 Jun 24;184(13):3573-3587.e29. doi: 10.1016/j.cell.2021.04.048. 
#PMID: 34062119; PMCID: PMC8238499. 

#Seurat Version 4.3.0.1
#SeuratObject Version 4.1.3
#R version 4.3.0 (2023-04-21) -- "Already Tomorrow"

#### Load Required Packages ####
library(Seurat)
library(RColorBrewer)
library(ggplot2)

#### Load Seurat Object ####
load("~/steele_et_al_2020_PDA.RData")

#### Integrate object by Sample ID--------------------------------------------------------------------------------------------------- ####
Idents(object = human_adj_norm_PDA) <- "ID"
human.list <- SplitObject(human_adj_norm_PDA, split.by = "ID")

human.list <- lapply(X = human.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

human.anchors <- FindIntegrationAnchors(object.list = human.list)
human_tissue_integrated <- IntegrateData(anchorset = human.anchors)
DefaultAssay(human_tissue_integrated) <- "integrated"

#Scale Data:
all.genes <- rownames(human_tissue_integrated)
human_tissue_integrated <- ScaleData(human_tissue_integrated, verbose = FALSE, features = all.genes)

#Run PCA and Determine Dimensions for 90% Variance:
human_tissue_integrated <- RunPCA(human_tissue_integrated, npcs = 30, verbose = FALSE)
st_dev <- human_tissue_integrated@reductions$pca@stdev
var <- st_dev^2
sum(var[1:14])/ sum(var)

#Find Neighbors and Cluster Cells:
human_tissue_integrated <- FindNeighbors(object = human_tissue_integrated, dims = 1:14)
human_tissue_integrated <- FindClusters(object = human_tissue_integrated, resolution = 2)

#Run UMAP:
human_tissue_integrated <- RunUMAP(human_tissue_integrated, reduction = "pca", dims = 1:14, verbose = F)
DimPlot(human_tissue_integrated, reduction = "umap", label = T)

#Identify Cell Populations:
DotPlot(human_tissue_integrated, features = c('KRT18', 'KRT8', 'CLU',
                                              'PRSS1', 'CTRB2',
                                              'MSLN', 'LRRN4',
                                              'COL1A2', 'PDPN',
                                              'CLEC3B', 'IL6',
                                              'TAGLN', 'ACTA2',
                                              'CD74',
                                              'CSPG4', 'RGS5',
                                              'CDH5', 'VWF', 'PLVAP',
                                              'POU2F3', 'DCLK1',
                                              "PTPRC", 
                                              "CD3E", "CD4", "FOXP3", "CD8A", "NKG7", 'PRF1',
                                              "CD68", "ARG1", "FCGR3A","ITGAM",'ITGAX',"CD14",'APOE',"MRC1", 
                                              "S100A8", "CD33",
                                              "HLA-DRB1","ITGAE","CLEC9A","BATF3", 'IRF7', 'CCL22',
                                              "CD79A", "CD19", "MS4A1","CCR10","PRDM1",
                                              "KIT","IL1RL1",
                                              "MKI67", 'CCNA2',
                                              'HBB'), cols = "RdBu", dot.scale = 8, assay = 'RNA') + RotatedAxis()

#Annotate Clusters:
Idents(human_tissue_integrated) <- "seurat_clusters"
human_tissue_integrated <- RenameIdents(human_tissue_integrated,
                                        "0" = "CD8+ T Lymphocyte", 
                                        "1" = "CD8+ T Lymphocyte", 
                                        "2" = "CD4+ T Lymphocyte", 
                                        "3" = "Macrophage", 
                                        "4" = "CK19+ ECAD+", 
                                        "5" = "CK19+ ECAD+", 
                                        "6" = "Fibroblast", 
                                        "7" = "Acinar", 
                                        "8" = "CK19+ ECAD+",
                                        "9" = "Granulocyte", 
                                        "10" = "CK19+ ECAD+", 
                                        "11" = "Macrophage", 
                                        "12" = "CK19+ ECAD+", 
                                        "13" = "Granulocyte", 
                                        "14" = "MasT Lymphocyte", 
                                        "15" = "CK19+ ECAD+", 
                                        "16" = "CK19+ ECAD+", 
                                        "17" = "Macrophage", 
                                        "18" = "CD4+ T Lymphocyte",
                                        "19" = "Granulocyte" , 
                                        "20" = "Naive T Lymphocyte",
                                        "21" = "Endocrine", 
                                        "22" = "Plasma Cell",
                                        "23" = "Granulocyte",
                                        "24" = "Macrophage", 
                                        "25" = "Pericyte",
                                        "26" = "Granulocyte", 
                                        "27" = "Treg", 
                                        "28" = "NK Cell", 
                                        "29" = "B Lymphocyte",
                                        "30" = "Dendritic Cell",
                                        "31" = "CK19+ ECAD+",
                                        "32" = "Dendritic Cell",
                                        "33" = "Acinar", 
                                        "34" = "Granulocyte",
                                        "35" = "RBC",
                                        "36" = "Granulocyte", 
                                        "37" = "CK19+ ECAD+ - Proliferating",
                                        "38" = "Acinar", 
                                        "39" = "Fibroblast", 
                                        "40" = "Endothelial", 
                                        "41" = "MasT Lymphocyte",
                                        "42" = "Treg",
                                        "43" = "CD8+ T Lymphocyte - Proliferating",
                                        "44" = "NK T Lymphocyte",
                                        "45" = "Macrophage - Proliferating", 
                                        "46" = "Macrophage",
                                        "47" = "Acinar",
                                        "48" = "Plasma Cell", 
                                        "49" = "Dendritic Cell",
                                        "50" = "Platelet", 
                                        "51" = "CK19+ ECAD+",
                                        "52" = "Plasma Cell", 
                                        "53" = "Granulocyte")
human_tissue_integrated[["manual_clusters"]] <- human_tissue_integrated@active.ident

Idents(human_tissue_integrated) <- "manual_clusters"
human_tissue_integrated <- RenameIdents(human_tissue_integrated,
                                        "CD8+ T Lymphocyte" = "CD8+ T Lymphocyte",  
                                        "CD4+ T Lymphocyte" = "CD4+ T Lymphocyte", 
                                        "Macrophage" = "Macrophage", 
                                        "CK19+ ECAD+" = "CK19+ ECAD+", 
                                        "Acinar" = "Acinar", 
                                        "Granulocyte" = "Granulocyte", 
                                        "MasT Lymphocyte" = "MasT Lymphocyte", 
                                        "Naive T Lymphocyte" = "Naive T Lymphocyte",
                                        "Endocrine" = "Endocrine", 
                                        "Plasma Cell" = "Plasma Cell", 
                                        "Pericyte" = "Pericyte",
                                        "Treg" = "CD4+ T Lymphocyte", 
                                        "NK Cell" = "NK Cell", 
                                        "B Lymphocyte" = "B Lymphocyte",
                                        "Dendritic Cell" = "Dendritic Cell",
                                        "RBC" = "RBC",
                                        "CK19+ ECAD+ - Proliferating" = "Proliferating",
                                        "Fibroblast" = "Fibroblast", 
                                        "Endothelial" = "Endothelial", 
                                        "CD8+ T Lymphocyte - Proliferating" = "Proliferating",
                                        "NK T Lymphocyte" = "NK T Lymphocyte",
                                        "Macrophage - Proliferating" = "Proliferating", 
                                        "Platelet" = "Platelet")

new_order <- c("CK19+ ECAD+", 
               "Fibroblast",
               "Pericyte",
               "Endothelial",
               "Macrophage",
               "Dendritic Cell",
               "Granulocyte",
               "CD4+ T Lymphocyte",
               "CD8+ T Lymphocyte",
               "Naive T Lymphocyte",
               "NK T Lymphocyte",
               "NK Cell",
               "B Lymphocyte",
               "Plasma Cell",
               "MasT Lymphocyte",
               "Acinar",
               "Endocrine",
               "Platelet",
               "RBC",
               "Proliferating")
human_tissue_integrated@active.ident <- factor(human_tissue_integrated@active.ident, levels = new_order)
human_tissue_integrated[["simple_clusters"]] <- human_tissue_integrated@active.ident

Idents(human_tissue_integrated) <- "simple_clusters"
human_tissue_integrated <- RenameIdents(human_tissue_integrated,
                                        "CD8+ T Lymphocyte" = "T Lymphocyte",  
                                        "CD4+ T Lymphocyte" = "T Lymphocyte", 
                                        "Macrophage" = "Myeloid", 
                                        "CK19+ ECAD+" = "CK19+ ECAD+", 
                                        "Acinar" = "Acinar", 
                                        "Granulocyte" = "Myeloid", 
                                        "MasT Lymphocyte" = "Myeloid", 
                                        "Naive T Lymphocyte" = "T Lymphocyte",
                                        "Endocrine" = "Endocrine", 
                                        "Plasma Cell" = "B Lymphocyte", 
                                        "Pericyte" = "Pericyte",
                                        "NK Cell" = "NK Cell", 
                                        "B Lymphocyte" = "B Lymphocyte",
                                        "Dendritic Cell" = "Myeloid",
                                        "RBC" = "RBC",
                                        "Proliferating" = "Proliferating",
                                        "Fibroblast" = "Fibroblast", 
                                        "Endothelial" = "Endothelial",
                                        "NK T Lymphocyte" = "T Lymphocyte",
                                        "Platelet" = "Platelet")

new_order <- c("CK19+ ECAD+", 
               "Fibroblast",
               "Pericyte",
               "Endothelial",
               "Myeloid",
               "T Lymphocyte",
               "NK Cell",
               "B Lymphocyte",
               "Acinar",
               "Endocrine",
               "Platelet",
               "RBC",
               "Proliferating")
human_tissue_integrated@active.ident <- factor(human_tissue_integrated@active.ident, levels = new_order)
human_tissue_integrated[["super_simple_clusters"]] <- human_tissue_integrated@active.ident

#Save Seurat Object:
save(human_tissue_integrated, file = 'huamn_integrated.RData')


#### Figures------------------------------------------------------------------------------------------------------------------------- ####
#### Figure 1B, Split Global UMAP #### 
Idents(human_tissue_integrated) <- "DiseaseState"
table(human_tissue_integrated@active.ident)
#PDAC AdjNorm 
#49018    9109 

Idents(human_tissue_integrated) <- "super_simple_clusters"
DimPlot(human_tissue_integrated, cols = c("CK19+ ECAD+" = "#fccde5", 
                                          "Fibroblast" = "#a6cee3",
                                          "Pericyte" = "#8dd3c7",
                                          "Endothelial" = "#1f78b4",
                                          "Myeloid" = "#cab2d6",
                                          "T Lymphocyte" = "#b2df8a",
                                          "NK Cell" = "#33a02c",
                                          "B Lymphocyte" = "#ffff99",
                                          "Acinar" = "#6a3d9a",
                                          "Endocrine" = "#fdbf6f",
                                          "Platelet" = "#ff7f00",
                                          "RBC" = "#b15928",
                                          "Proliferating" = "grey"), split.by = "DiseaseState")

#### Figure S1A, Global Cluster Marker Dotplot #### 
Idents(human_tissue_integrated) <- "super_simple_clusters"
DotPlot(human_tissue_integrated, features = c('KRT18', 'KRT8', "KRT19", "CDH1",
                                              'COL1A2', 'PDPN',
                                              'CSPG4', 'RGS5',
                                              'CDH5', 'VWF', 'PLVAP',
                                              "PTPRC",
                                              "CD68", "FCGR3A","CD14","HLA-DRB1", 
                                              "CD3E", "CD4", "CD8A",
                                              "NKG7", 'PRF1',
                                              "CD79A", "CD19", "MS4A1","CCR10",
                                              "PRSS1", "CTRB2",
                                              "INS", "SST", 'GCG',
                                              "ITGA2B", "ITGB3", "SELP",
                                              'HBA1', "HBB",
                                              "MKI67", 'CCNA2'), cols = "RdYlBu", dot.scale = 8, assay = 'RNA') + RotatedAxis()

#### Figure 1C, IL33 FeaturePlots ####
Idents(human_tissue_integrated) <- "DiseaseState"
tissue_PDA <- subset(human_tissue_integrated, idents = "PDAC")
tissue_norm <- subset(human_tissue_integrated, idents = "AdjNorm")
DefaultAssay(tissue_PDA) <- "RNA"
DefaultAssay(tissue_norm) <- "RNA"

FeaturePlot(tissue_PDA, features = "IL33", order = T) + scale_colour_gradientn(colours = brewer.pal(n = 9, name = "Reds"))
FeaturePlot(tissue_norm, features = "IL33", order = T) + scale_colour_gradientn(colours = brewer.pal(n = 9, name = "Reds"))