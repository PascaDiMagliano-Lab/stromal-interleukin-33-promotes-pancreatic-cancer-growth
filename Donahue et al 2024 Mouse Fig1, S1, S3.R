#This script recreates the murine scRNAseq analysis from Figures 1D&E, S1B&C, and S3C of Donahue et al. 2024.
#This analysis utilizes datasets generated in this study, as well as previously published datasets from:

#Healthy 1&2 (NCBI GEO GSM5011580 & GSM5011581)
#Kemp, S.B., et al., Apolipoprotein E Promotes Immune Suppression in Pancreatic Cancer through NF-kappaB-Mediated Production of CXCL1. 
#Cancer Res, 2021. 81(16): p. 4305-4318. 

#Healthy 3 (NCBI GEO GSM3577882)
#Hosein, A.N., et al., Cellular heterogeneity during mouse pancreatic ductal adenocarcinoma progression at single-cell resolution. 
#JCI Insight, 2019. 5(16).

#PDA 1 (NCBI GEO GSM6127792)
#Du, W., et al., WNT signaling in the tumor microenvironment promotes immunosuppression in murine pancreatic cancer.
#J Exp Med, 2023. 220(1).

#PDA 2 (NCBI GEO GSE129455)
#Elyada, E., et al., Cross-Species Single-Cell Analysis of Pancreatic Ductal Adenocarcinoma Reveals Antigen-Presenting Cancer-Associated Fibroblasts. 
#Cancer Discov, 2019. 9(8): p. 1102-1123.

#PDA 3 (NCBI GEO GSM3577885)
#Hosein, A.N., et al., Cellular heterogeneity during mouse pancreatic ductal adenocarcinoma progression at single-cell resolution. 
#JCI Insight, 2019. 5(16).

#PanIN 1-3 (this manuscript, NCBI GEO GSE269182 and GSE269888)

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
library(dplyr)
library(stringr)
library(RColorBrewer)
library(ggplot2)
library(scCustomize)
library(msigdbr)
library(fgsea)

#### Preprocessing ------------------------------------------------------------------------------------------------------------ ####

#Load Raw Data:

# HEALTHY #
healthy_data1 <- Read10X("~/Kemp_et_al_GSM5011580/raw_feature_bc_matrix/")
healthy_data2 <- Read10X_h5("~/Kemp_et_al_GSM5011581/raw_feature_bc_matrix.h5")
healthy_data3 <- Read10X_h5("~/Hosein_et_al_GSM3577882/filtered_gene_bc_matrices_h5.h5")

# PANIN #
KC_data_769 <- Read10X_h5("~/KC 769/raw_feature_bc_matrix.h5")
KC_data_4751 <- Read10X_h5("~/KC 4751/raw_feature_bc_matrix.h5")
KC_data_5088 <- Read10X_h5("~/KC 5088/raw_feature_bc_matrix.h5")

# PDA #
PDA_data1 <- Read10X_h5("~/Du_et_al_GSM6127792/filtered_feature_bc_matrix.h5")
PDA_data2 <- Read10X_h5("~/Elyada_et_al_GSE129455/filtered_feature_bc_matrix.h5")
PDA_data3 <- Read10X_h5("~/Hosein_et_al_GSM3577885/filtered_gene_bc_matrices_h5.h5")

#Create Seurat Objects:
healthy1 <- CreateSeuratObject(healthy_data1, min.cells = 3, min.features = 100)
healthy2 <- CreateSeuratObject(healthy_data2, min.cells = 3, min.features = 100)
healthy3 <- CreateSeuratObject(healthy_data3, min.cells = 3, min.features = 100)

KC_769 <- CreateSeuratObject(KC_data_769, min.cells = 3, min.features = 100)
KC_4751 <- CreateSeuratObject(KC_data_4751, min.cells = 3, min.features = 100)
KC_5088 <- CreateSeuratObject(KC_data_5088, min.cells = 3, min.features = 100)

PDA1 <- CreateSeuratObject(PDA_data1, min.cells = 3, min.features = 100)
PDA2 <- CreateSeuratObject(PDA_data2, min.cells = 3, min.features = 100)
PDA3 <- CreateSeuratObject(PDA_data3, min.cells = 3, min.features = 100)

#Add Desired Metadata:

healthy1[["Disease"]] <- "Healthy"
healthy2[["Disease"]] <- "Healthy"
healthy3[["Disease"]] <- "Healthy"

KC_769[["Disease"]] <- "PanIN"
KC_4751[["Disease"]] <- "PanIN"
KC_5088[["Disease"]] <- "PanIN"

PDA1[["Disease"]] <- "PDA"
PDA2[["Disease"]] <- "PDA"
PDA3[["Disease"]] <- "PDA"

healthy1[["Sample"]] <- "Healthy_1"
healthy2[["Sample"]] <- "Healthy_2"
healthy3[["Sample"]] <- "Healthy_3"

KC_769[["Sample"]] <- "PanIN_1"
KC_4751[["Sample"]] <- "PanIN_2"
KC_5088[["Sample"]] <- "PanIN_3"

PDA1[["Sample"]] <- "PDA_1"
PDA2[["Sample"]] <- "PDA_2"
PDA3[["Sample"]] <- "PDA_3"

#Merge Seurat Objects:
murine_collab <- merge(x = healthy1, y = c(healthy2, 
                                           healthy3, 
                                           KC_769,
                                           KC_4751,
                                           KC_5088,
                                           PDA1,
                                           PDA2,
                                           PDA3), add.cell.ids = (c("a","b","c","d", "e", "f", "g", "h", "i")))

#Normalize Data:
murine_collab <- NormalizeData(object = murine_collab, normalization.method = "LogNormalize", scale.factor = 10000)

#Apply Unbiased QC Cutoffs:
murine_collab[["percent.mt"]] <- PercentageFeatureSet(object = murine_collab, pattern = "^mt-")
murine_collab <- subset(x = murine_collab, subset = nCount_RNA > 1000 & nCount_RNA < 100000 & percent.mt < 15)

#Integrate object by Sample:
Idents(object = murine_collab) <- "Sample"
murine.list <- SplitObject(murine_collab, split.by = "Sample")

murine.list <- lapply(X = murine.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

murine.anchors <- FindIntegrationAnchors(object.list = murine.list)
murine_collab <- IntegrateData(anchorset = murine.anchors)
DefaultAssay(murine_collab) <- "integrated"

#Scale Data:
all.genes <- rownames(murine_collab)
murine_collab <- ScaleData(murine_collab, verbose = FALSE, features = all.genes)

#Run PCA and Determine Dimensions for 90% Variance:
murine_collab <- RunPCA(murine_collab, npcs = 30, verbose = FALSE)
st_dev <- murine_collab@reductions$pca@stdev
var <- st_dev^2
sum(var[1:16])/ sum(var)

#Find Neighbors and Cluster Cells:
murine_collab <- FindNeighbors(object = murine_collab, dims = 1:16)
murine_collab <- FindClusters(object = murine_collab, resolution = 2.8)

#Run UMAP:
murine_collab <- RunUMAP(murine_collab, reduction = "pca", dims = 1:16, verbose = F)
DimPlot(murine_collab, reduction = "umap", label = T)
DimPlot(murine_collab, reduction = "umap", label = T, group.by = "seurat_clusters")

#Identify Cell Populations:
Idents(murine_collab) <- "seurat_clusters"
DotPlot(murine_collab, features = c("Krt19", "Cdh1",
                                    "Try4", "Amy2a2", 
                                    "Msln", 
                                    "Col1a2", "Pdpn", 
                                    "Clec3b", "Il6",
                                    "Tagln", "Acta2", 
                                    "Saa3", "Cd74",  
                                    "Cspg4", 
                                    "Pecam1", "Cdh5",  
                                    "Pou2f3", "Dclk1", 
                                    "Ptprc",  
                                    "Cd68", "Adgre1", "Itgam", "Cd14", "Mrc1", 
                                    "S100a8", "Cd33", "Ly6c2","Ly6g", 
                                    "H2-Eb1", "Itgae", "Clec9a", "Batf3", 
                                    "Cd3e", "Cd4", "Foxp3","Cd8a", 
                                    "Trdc", 
                                    "Nkg7", "Klrg1", 
                                    "Il1rl1", "Gata3", 
                                    "Cd79a", 
                                    "Cd19", "Ms4a1", 
                                    "Ccr10", "Prdm1", 
                                    "Kit",
                                    "Mki67", "Ccna2", 
                                    "Hbb-bt"), cols = "RdYlBu", dot.scale = 8, assay = 'RNA') + RotatedAxis()

DotPlot(murine_collab, features = c("Krt18","Krt8", "Tff1", "Lyz", 
                                    "Mmp7", "Clu", "Spp1", 
                                    "Krt19", "Spink1",
                                    "Cdh1", "Epcam" ), cols = "RdBu", dot.scale = 5, assay = 'RNA') + RotatedAxis()

#Annotate Clusters:
Idents(murine_collab) <- "seurat_clusters"
murine_collab <- RenameIdents(murine_collab,
                              "0" = "Fibroblast", 
                              "1" = "Macrophage", 
                              "2" = "Granulocyte", 
                              "3" = "Fibroblast", 
                              "4" = "Macrophage", 
                              "5" = "Fibroblast", 
                              "6" = "Fibroblast", 
                              "7" = "Macrophage", 
                              "8" = "Fibroblast",
                              "9" = "B Cell", 
                              "10" = "Macrophage", 
                              "11" = "Granulocyte", 
                              "12" = "T Cell", 
                              "13" = "Macrophage", 
                              "14" = "Fibroblast", 
                              "15" = "Macrophage", 
                              "16" = "Dendritic Cell", 
                              "17" = "Proliferating Macrophage",
                              "18" = "Fibroblast", 
                              "19" = "T Cell", 
                              "20" = "Macrophage", 
                              "21" = "CK19+ ECAD+",
                              "22" = "Acinar", 
                              "23" = "RBC", 
                              "24" = "Fibroblast", 
                              "25" = "Granulocyte", 
                              "26" = "Acinar", 
                              "27" = "Acinar", 
                              "28" = "Fibroblast",
                              "29" = "Plasma",
                              "30" = "CK19+ ECAD+",
                              "31" = "CK19+ ECAD+",
                              "32" = "CK19+ ECAD+",
                              "33" = "T Cell",
                              "34" = "Macrophage", 
                              "35" = "Endothelial", 
                              "36" = "Macrophage",
                              "37" = "Fibroblast",
                              "38" = "T Cell",
                              "39" = "Fibroblast",
                              "40" = "RBC",
                              "41" = "Mesothelial",
                              "42" = "Acinar", 
                              "43" = "B Cell",
                              "44" = "Dendritic Cell",
                              "45" = "Mesothelial", 
                              "46" = "Tuft Cell", 
                              "47" = "CK19+ ECAD+", 
                              "48" = "RBC", 
                              "49" = "Fibroblast", 
                              "50" = "Proliferating Fibroblast", 
                              "51" = "Proliferating T Cell",
                              "52" = "Acinar",
                              "53" = "Monocyte",
                              "54" = "Mast Cell",
                              "55" = "B Cell")
murine_collab[["manual_clusters"]] <- murine_collab@active.ident

Idents(murine_collab) <- "manual_clusters"
murine_collab <- RenameIdents(murine_collab,
                              "CK19+ ECAD+" = "CK19+ ECAD+", 
                              "Acinar" = "Acinar",
                              "Tuft Cell" = "Tuft Cell",
                              "Endothelial" = "Endothelial",
                              "Mesothelial" = "Mesothelial",
                              "Fibroblast" = "Fibroblast",
                              "Macrophage" = "Myeloid",
                              "Monocyte" = "Myeloid",
                              "Dendritic Cell" = "Myeloid",
                              "Granulocyte" = "Myeloid",
                              "Mast Cell" = "Myeloid",
                              "T Cell" = "Lymphocyte",
                              "B Cell" = "Lymphocyte",
                              "Plasma" = "Lymphocyte",
                              "RBC" = "RBC",
                              "Proliferating Macrophage" = "Proliferating",
                              "Proliferating Fibroblast" = "Proliferating",
                              "Proliferating T Cell" = "Proliferating")
new_order <- c("CK19+ ECAD+", 
               "Acinar",
               "Tuft Cell",
               "Endothelial",
               "Mesothelial",
               "Fibroblast",
               "Myeloid",
               "Lymphocyte",
               "RBC",
               "Proliferating")
murine_collab@active.ident <- factor(murine_collab@active.ident, levels = new_order)
murine_collab[["simple_clusters"]] <- murine_collab@active.ident

#Save Seurat Object:
save(murine_collab, file = 'UMich_and_collab_mouse_atlas.RData')

#### Figures------------------------------------------------------------------------------------------------------------------------- ####

plot_colors <- c("CK19+ ECAD+" = "gold", 
                 "Acinar" = "tan1",
                 "Tuft Cell" = "violetred4",
                 "Endothelial" = "#6a3d9a",
                 "Mesothelial" = "darkblue",
                 "Fibroblast" = "skyblue2",
                 "Myeloid" = "#7fbc41",
                 "Lymphocyte" = "#fec5e5",
                 "NK Cell" = "hotpink2",
                 "RBC" = "snow3",
                 "Proliferating" = "honeydew2")

#Version of Object with Proliferating and Red Blood Cell Populations Excluded:
murine_collab_plots <- subset(murine_collab, idents = c("Proliferating", "RBC"), invert = T)

#### Figure 1D, Split Global UMAP #### 
DimPlot(murine_collab, reduction = "umap", label = F, 
        cols = plot_colors, split.by = "Disease", pt.size = 0.4)

Idents(murine_collab_plots) <- "Disease"
table(murine_collab_plots@active.ident)
#Healthy   PanIN     PDA 
#   6336   14207   18834

#### Figure S1B, Global Cluster Marker Dotplot ####
DotPlot(murine_collab, features = c("Krt19", "Cdh1",
                                    "Try4", "Amy2a2",
                                    "Pou2f3", "Dclk1",
                                    "Pecam1", "Cdh5",
                                    "Msln", "Lrrn4",
                                    "Col1a2", "Pdpn", "Pdgfra", "Pdgfrb",
                                    "Ptprc",
                                    "Cd68", "Adgre1", "Itgam", "Cd14",
                                    "H2-Eb1", "Batf3",
                                    "S100a8", 
                                    "Cd3e", "Cd4", "Cd8a",
                                    "Cd79a", "Cd19", "Ms4a1", "Ccr10",
                                    "Nkg7", "Gzmb",
                                    "Hbb-bs", "Hba-a1", 
                                    "Mki67", "Ccna2"), dot.scale = 8, cols = 'BrBG', assay = 'RNA') + RotatedAxis()

#### Figure 1E, Il33 Dotplot #### 
DotPlot(murine_collab_plots, features = "Il33", dot.scale = 8, cols = 'BrBG', assay = 'RNA', split.by = "Disease", group.by = "simple_clusters") + RotatedAxis()

#### Figure S1C, Cell Abundance Bar Graph ####
Idents(object = murine_collab) <- 'Sample'
table(murine_collab@active.ident)
#Healthy_1   Healthy_2   Healthy_3   PanIN_1   PanIN_2   PanIN_3  PDA_1  PDA_2 PDA_3 
#     2590        1787        2482      5325      4448      5368   4999  14029   872 

plot_IDs <- levels(murine_collab)
samples.list <- unique(murine_collab$Sample)
clusters <- lapply(samples.list, function(x){
  subset <- subset(murine_collab, subset = Sample == x)
  dist <- data.frame(table(subset$simple_clusters))
  
  return(dist)
})

names(clusters) <- samples.list

clusters_percent <- lapply(clusters, FUN = function(x){
  summ <- sum(x$Freq)
  x$Freq <- (x$Freq/summ)
  return(x)
})

clusters_dist <- reshape2::melt(clusters, id.var = "Var1")
colnames(clusters_dist) <- c("simple_clusters","variable","value","Sample")
clusters_percent_dist <- reshape2::melt(clusters_percent, id.var = "Var1")
colnames(clusters_percent_dist) <- c("simple_clusters","variable","value","Sample")

ggplot(clusters_percent_dist, aes(fill=simple_clusters, y = value, x = Sample)) + 
  scale_x_discrete(limits = plot_IDs) + 
  geom_bar(position = "stack", stat = "identity") + 
  theme_bw()+
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, size = 8, face = "bold", hjust = 1))+
  ylab("count") + xlab("orig.ident") + ggtitle("Relative Cell Types Abundance") +scale_fill_manual(values = plot_colors)

#### Figure S3C, Hallmark ROS Violin Plots ####
h_gene_sets = msigdbr(species = "mouse", category = "H")
msigdbr_list = split(x = h_gene_sets$gene_symbol, f = h_gene_sets$gs_name)

ros <- msigdbr_list[["HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY"]]

DefaultAssay(murine_collab_plots) <- "RNA"
murine_collab_plots <- AddModuleScore(murine_collab_plots,
                                      features = list(ros),
                                      name="Hallmark_ROS_Pathway")

#Plot module scores:
VlnPlot(murine_collab_plots, features = "Hallmark_ROS_Pathway1", split.by = "Disease", group.by = "simple_clusters", pt.size = 0, cols = c("lavender", "#c097c5", "#69449b"))+ geom_hline(yintercept=0, linetype="dashed", 
                                                                                                                                                                                          linewidth=0.5)
