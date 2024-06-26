#This script recreates the CreER and CreER;Il33f/f orthotopic scRNAseq analysis from 
#Figures 2H-J, 3A-H, 4B-C, S2D-H, S3D, S4A, S4G, S5B-D of Donahue et al. 2024.

#Data is available at NCBI GEO GSE269888.

#Data was processed in line with the Seurat workflow:
#Website: https://satijalab.org/seurat/index.html

#Reference: Hao et al., Integrated analysis of multimodal single-cell data. 
#Cell. 2021 Jun 24;184(13):3573-3587.e29. doi: 10.1016/j.cell.2021.04.048. 
#PMID: 34062119; PMCID: PMC8238499. 

#Cell Interaction Analysis performed in line with the CellChat workflow (version 1.6.1)
#Website: https://github.com/sqjin/CellChat

#Reference: Jin, S., et al., Inference and analysis of cell-cell communication using CellChat. 
#Nature Communications, 2021. 12(1): p. 1088.

#Seurat Version 4.3.0.1
#SeuratObject Version 4.1.3
#CellChat Version 1.6.1
#R version 4.3.0 (2023-04-21) -- "Already Tomorrow"

#### Load Required Packages ####
library(Seurat)
library(dplyr)
library(stringr)
library(RColorBrewer)
library(ggplot2)
library(scCustomize)
library(viridis)
library(fgsea)
library(msigdbr)
library(EnhancedVolcano)
library(ComplexHeatmap)
library(CellChat)
library(patchwork)

#### Preprocessing - Global Object ----------------------------------------------------------------------------------------------------------- ####

#Load Raw Data:
wt_data <- Read10X_h5("~/Dropbox (University of Michigan)/2. R/2022 Katelyn's IL33 paper/OT3/WT 7940b OT/filtered_feature_bc_matrix.h5")
ko_data <- Read10X_h5("~/Dropbox (University of Michigan)/2. R/2022 Katelyn's IL33 paper/OT3/7940b PDGFRa;IL33 OT/filtered_feature_bc_matrix.h5")

#Create Seurat Objects:
wt <- CreateSeuratObject(wt_data, min.cells = 3, min.features = 100)
ko <- CreateSeuratObject(ko_data, min.cells = 3, min.features = 100)

#Add Desired Metadata:

wt[["Group"]] <- "Pdgfra-CreERT2"
ko[["Group"]] <- "Pdgfra-CreERT2;IL33 f/f"

wt[["Run_Date"]] <- "Jan_2022"
ko[["Run_Date"]] <- "Jan_2022"

#Merge Seurat Objects:
IL33_OT <- merge(x = wt, y = ko, add.cell.ids = (c("wt","ko")))

#Normalize Data:
IL33_OT <- NormalizeData(object = IL33_OT, normalization.method = "LogNormalize", scale.factor = 10000)

#Apply Unbiased QC Cutoffs:
IL33_OT[["percent.mt"]] <- PercentageFeatureSet(object = IL33_OT, pattern = "^mt-")
IL33_OT <- subset(x = IL33_OT, subset = nCount_RNA > 1000 & nCount_RNA < 100000 & percent.mt < 15)

#Find Variable Genes:
IL33_OT <- FindVariableFeatures(IL33_OT, selection.method = "vst", nfeatures = 2000)

#Scale Data:
all.genes <- rownames(IL33_OT)
IL33_OT <- ScaleData(IL33_OT, verbose = T, features = all.genes)

#Run PCA and Determine Dimensions for 90% Variance:
IL33_OT <- RunPCA(IL33_OT, npcs = 30, verbose = FALSE)
st_dev <- IL33_OT@reductions$pca@stdev
var <- st_dev^2
sum(var[1:19])/ sum(var)

#Find Neighbors and Cluster Cells:
IL33_OT <- FindNeighbors(object = IL33_OT, dims = 1:19)
IL33_OT <- FindClusters(object = IL33_OT, resolution = 7)

#Run UMAP:
IL33_OT <- RunUMAP(IL33_OT, reduction = "pca", dims = 1:19, verbose = F)
DimPlot(IL33_OT, reduction = "umap", label = T)

#Identify Cell Populations:
DotPlot(IL33_OT, features = c("Krt19", "Epcam","Msln","Try4", "Col1a2", "Acta2", "Clec3b", "Cspg4","Saa3","Ptprc", "Cd3e", "Trdc","Cd4", "Foxp3", "Ccr7", "Cd8a", "Nkg7","Il2ra","Il1rl1","Arg1", "Ccr2","Cd68","Fcgr3","Adgre1", "Itgam","Cd14","Mrc1", "S100a8", "Cd33","H2-Eb1",
                                        "Itgae","Clec9a","Batf3", "Ly6c2", "Ly6g","Cd79a", "Cd19", "Ms4a1","Ccr10","Prdm1","Kit","Ccna2", "Mki67","Pecam1", "Cdh5", "Hbb-bt"), cols = "RdYlBu", dot.scale = 8) + RotatedAxis()

DotPlot(IL33_OT, features = c("Krt19", "Pdpn", "Pdgfra","Ero1l","Vegfa"), cols = "RdYlBu", dot.scale = 8) + RotatedAxis() #EMT Cells
DotPlot(IL33_OT, features = c("Krt19", "Epcam","Cdh1","Msln","Lrrn4","Pdgfra", "Col1a2"), cols = "RdYlBu", dot.scale = 8) + RotatedAxis() #Mesothelial
DotPlot(IL33_OT, features = c("Mki67", "Ccna2"), cols = "RdYlBu", dot.scale = 8) + RotatedAxis() #Proliferating

#Annotate Clusters:
Idents(IL33_OT) <- "seurat_clusters"
IL33_OT <- RenameIdents(IL33_OT,
                                  "0" = "CD4+ T Cell", 
                                  "1" = "Tumor", 
                                  "2" = "MonoMac", 
                                  "3" = "Proliferating", 
                                  "4" = "MonoMac", 
                                  "5" = "CD8+ T Cell", 
                                  "6" = "Granulocyte", 
                                  "7" = "CD4+ T Cell", 
                                  "8" = "CD4+ T Cell",
                                  "9" = "MonoMac", 
                                  "10" = "MonoMac", 
                                  "11" = "Tumor", 
                                  "12" = "Fibroblast", 
                                  "13" = "Macrophage", 
                                  "14" = "Macrophage", 
                                  "15" = "Macrophage", 
                                  "16" = "Granulocyte", 
                                  "17" = "Macrophage",
                                  "18" = "Tumor", 
                                  "19" = "Fibroblast", 
                                  "20" = "Proliferating", 
                                  "21" = "Tumor",
                                  "22" = "B Cell", 
                                  "23" = "Dendritic Cell", 
                                  "24" = "Fibroblast", 
                                  "25" = "Granulocyte", 
                                  "26" = "Tumor", 
                                  "27" = "Macrophage", 
                                  "28" = "Granulocyte",
                                  "29" = "Macrophage",
                                  "30" = "Macrophage",
                                  "31" = "Tumor",
                                  "32" = "CD8+ T Cell",
                                  "33" = "Tumor",
                                  "34" = "Plasma Cell", 
                                  "35" = "MonoMac", 
                                  "36" = "Fibroblast",
                                  "37" = "Acinar",
                                  "38" = "Fibroblast",
                                  "39" = "B Cell",
                                  "40" = "Macrophage",
                                  "41" = "MonoMac",
                                  "42" = "Mesothelial", 
                                  "43" = "Tumor",
                                  "44" = "Fibroblast",
                                  "45" = "Endothelial",
                                  "46" = "Tumor",
                                  "47" = "Dendritic Cell",
                                  "48" = "Macrophage", 
                                  "49" = "Dendritic Cell", 
                                  "50" = "Macrophage",
                                  "51" = "Proliferating",
                         "52" = "Plasma Cell",
                         "53" = "Proliferating",
                         "54" = "Acinar",
                         "55" = "RBC",
                         "56" = "Dendritic Cell",
                         "57" = "MonoMac", 
                         "58" = "Fibroblast",
                         "59" = "Proliferating",
                         "60" = "RBC",
                         "61" = "Granulocyte",
                         "62" = "Granulocyte",
                         "63" = "Macrophage", 
                         "64" = "RBC", 
                         "65" = "Proliferating",
                         "66" = "Macrophage",
                         "67" = "RBC",
                         "68" = "Ductal",
                         "69" = "Macrophage",
                         "70" = "RBC",
                         "71" = "Pericyte",
                         "72" = "Tumor", 
                         "73" = "Proliferating",
                         "74" = "Endothelial",
                         "75" = "Proliferating",
                         "76" = "Proliferating",
                         "77" = "Proliferating",
                        "78" = "CD4+ T Cell")

#Manually Identify NK Cells and ILC2s:
t_cells <- subset(IL33_OT, idents = c("CD4+ T Cell", "CD8+ T Cell"))
nk_idents <- WhichCells(t_cells, expression = Gzmb > 0 & Nkg7 > 0 & Cd3e == 0 & Cd8a == 0 & Trdc == 0)
ilc2_idents <- WhichCells(t_cells, expression = Gata3 > 0 & Il1rl1 > 0 & Cd3e == 0 & Cd8a == 0 & Cd4 == 0 & Foxp3 == 0)

#Manually Identify Mast Cells:
myeloid <- subset(IL33_OT, idents = c("MonoMac", "Macrophage", "Dendritic Cell"))
mast_idents <- WhichCells(myeloid, expression = Il1rl1 > 0 & Kit > 0 & Adgre1 == 0 & Batf3 == 0)

#Manually Identify EMT-Like Cells:
emt_meso <- subset(IL33_OT, idents = c("Fibroblast", "Mesothelial"))
emt_idents <- WhichCells(emt_meso, expression = Krt19 > 0 & Ero1l > 1 & Vegfa > 1 & Lrrn4 == 0 & Msln == 0)

#Apply Manual Annotations:
Idents(object = IL33_OT, cells = nk_idents) <- "NK Cell"
Idents(object = IL33_OT, cells = ilc2_idents) <- "ILC2"
Idents(object = IL33_OT, cells = mast_idents) <- "Mast Cell"
Idents(object = IL33_OT, cells = emt_idents) <- "EMT"

#Organize Identities:
new_order <- c("Tumor",
               "EMT",
               "Ductal",
               "Acinar",
               "Endothelial",
               "Mesothelial",
               "Fibroblast",
               "Pericyte",
               "Macrophage",
               "MonoMac",
               "Dendritic Cell",
               "Granulocyte",
               "CD8+ T Cell",
               "NK Cell",
               "CD4+ T Cell",
               "ILC2",
               "Mast Cell",
               "B Cell",
               "Plasma Cell",
               "Proliferating",
               "RBC")
IL33_OT@active.ident <- factor(IL33_OT@active.ident, levels = new_order)
IL33_OT[["paper_clusters"]] <- IL33_OT@active.ident

#Save Seurat Object:
save(IL33_OT, file = 'IL33_OT.RData')

#### Preprocessing - ST2+ Cells -------------------------------------------------------------------------------------------------------------- ####
#Subset and Re-Cluster ST2-Potential Cell Populations Alone:
st2 <- subset(IL33_OT, idents = c("CD8+ T Cell", "CD4+ T Cell", "NK Cell","Mast Cell", "ILC2"))

#Find Variable Genes:
st2 <- FindVariableFeatures(st2, selection.method = "vst", nfeatures = 2000)

#Scale Data:
all.genes <- rownames(st2)
st2 <- ScaleData(st2, verbose = T, features = all.genes)

#Run PCA and Determine Dimensions for 90% Variance:
st2 <- RunPCA(st2, npcs = 30, verbose = FALSE)
st_dev <- st2@reductions$pca@stdev
var <- st_dev^2
sum(var[1:23])/ sum(var)

#Find Neighbors and Cluster Cells:
st2 <- FindNeighbors(object = st2, dims = 1:23)
st2 <- FindClusters(object = st2, resolution = 4.5)

#Run UMAP:
st2 <- RunUMAP(st2, reduction = "pca", dims = 1:23, verbose = F)
DimPlot(st2, reduction = "umap", label = T)

#Identify Cell Populations:
DotPlot(st2, 
        features = c("Krt19", "Krt8","Krt18","Epcam",
                     "Ero1l", "Vegfa", "Twist1", "Snai1",
                     "Mmp7", "Clu", "Spp1", 
                     "Try4", "Amy2a2", 
                     "Pecam1", "Cdh5",
                     "Msln", "Lrrn4",
                     "Col1a1","Pdpn", "Pdgfra", "Pdgfrb",
                     "Cspg4",
                     "Ptprc",
                     "Cd68","Itgam", "Cd14", "Adgre1", "Mrc1",
                     "Ccr2","Csf1r", "Tlr2",
                     "H2-Eb1","Itgae","Clec9a","Batf3",
                     "S100a8", "Cd33", 
                     "Cd3e", "Cd8a","Trdc",
                     "Nkg7", "Gzmb",
                     "Cd4", 
                     "Gata3", "Il1rl1",
                     "Mitf", "Kit",
                     "Cd79a", "Cd19", "Ms4a1","Ccr10","Prdm1",
                     "Mki67", "Ccna2",
                     "Hbb-bt", "Hbb-bs"), cols = "RdYlBu", dot.scale = 8) + RotatedAxis()

#Remove CD8+ T cells, NK cells, and Contaminating Macrophages:
st2 <- subset(st2, idents = c(0, 1, 2, 3, 8, 18, 20), invert = T)

#Find Variable Genes:
st2 <- FindVariableFeatures(st2, selection.method = "vst", nfeatures = 2000)

#Scale Data:
all.genes <- rownames(st2)
st2 <- ScaleData(st2, verbose = T, features = all.genes)

#Run PCA and Determine Dimensions for 90% Variance:
st2 <- RunPCA(st2, npcs = 30, verbose = FALSE)
st_dev <- st2@reductions$pca@stdev
var <- st_dev^2
sum(var[1:24])/ sum(var)

#Find Neighbors and Cluster Cells:
st2 <- FindNeighbors(object = st2, dims = 1:24)
st2 <- FindClusters(object = st2, resolution = 4.5)

#Run UMAP:
st2 <- RunUMAP(st2, reduction = "pca", dims = 1:24, verbose = F)
DimPlot(st2, reduction = "umap", label = T)

#Remove Gamma Delta T Cells and Remaining CD8+ T cells:
st2 <- subset(st2, idents = c(0, 1, 6, 15), invert = T)

#Find Variable Genes:
st2 <- FindVariableFeatures(st2, selection.method = "vst", nfeatures = 2000)

#Scale Data:
all.genes <- rownames(st2)
st2 <- ScaleData(st2, verbose = T, features = all.genes)

#Run PCA and Determine Dimensions for 90% Variance:
st2 <- RunPCA(st2, npcs = 30, verbose = FALSE)
st_dev <- st2@reductions$pca@stdev
var <- st_dev^2
sum(var[1:25])/ sum(var)

#Find Neighbors and Cluster Cells:
st2 <- FindNeighbors(object = st2, dims = 1:25)
st2 <- FindClusters(object = st2, resolution = 2.5)

#Run UMAP:
st2 <- RunUMAP(st2, reduction = "pca", dims = 1:25, verbose = F)
DimPlot(st2, reduction = "umap", label = T)

#Remove Remaining CD8+ T cells:
st2 <- subset(st2, idents = 4, invert = T) 

#Find Variable Genes:
st2 <- FindVariableFeatures(st2, selection.method = "vst", nfeatures = 2000)

#Scale Data:
all.genes <- rownames(st2)
st2 <- ScaleData(st2, verbose = T, features = all.genes)

#Run PCA and Determine Dimensions for 90% Variance:
st2 <- RunPCA(st2, npcs = 30, verbose = FALSE)
st_dev <- st2@reductions$pca@stdev
var <- st_dev^2
sum(var[1:25])/ sum(var)

#Find Neighbors and Cluster Cells:
st2 <- FindNeighbors(object = st2, dims = 1:25)
st2 <- FindClusters(object = st2, resolution = 2.8)

#Run UMAP:
st2 <- RunUMAP(st2, reduction = "pca", dims = 1:25, verbose = F)
DimPlot(st2, reduction = "umap", label = T)

#Identify Cell Populations:
DotPlot(st2, 
        features = c("Il1rl1", "Cd4","Cd8a","Cd3e","Foxp3", "Gata3", "Il17a", "Tbx21", "Cma1", "Arg1", "Il4", "Il13"), cols = "RdYlBu", dot.scale = 8) + RotatedAxis()

DotPlot(st2, features = c("Cd3e", "Cd3d","Sell", "Ccr7", "Tcf7","Il7r", "Cd28", "Cd4", "Tbx21", "Ifng", "Cxcr3", "Arg1","Gata3", "Il1rl1", "Il4", "Rorc", "Il17a", "Il22", "Il2ra", "Foxp3", "Icos", "Cd8a", "Gzmb","Prf1", "Eomes","Lag3","Pdcd1","Trdc","Nkg7", "Klrg1", "Ccna2", "Mki67"), cols = "RdYlBu")+ RotatedAxis()

#Annotate Clusters:
Idents(st2) <- "seurat_clusters"
st2 <- RenameIdents(st2,
                    "0" = "Treg", 
                    "1" = "Treg", 
                    "2" = "Th1", 
                    "3" = "Naive CD4", 
                    "4" = "ILC2", 
                    "5" = "Treg", 
                    "6" = "Treg", 
                    "7" = "Th17", 
                    "8" = "Exhausted CD4",
                    "9" = "Exhausted CD4", 
                    "10" = "Exhausted CD4",
                    "11" = "Th1")

#Manually Identify Mast Cells:
mast_ids <- WhichCells(st2, idents = c("ILC2", "Th1"), expression = Kit > 0 & Cma1 > 0 | Kit > 0 & Mitf > 0 )
FeaturePlot(st2, features = "Mitf")
FeaturePlot(st2, features = "Kit")

#Apply Manual Annotations:
Idents(st2, cells = mast_ids) <- "Mast Cell"

#Save Annotations:
new_order <- c("Naive CD4",
               "Exhausted CD4",
               "Th1",
               "Th17",
               "Treg",
               "ILC2",
               "Mast Cell")
st2@active.ident <- factor(st2@active.ident, levels = new_order)
st2[["st2_clusters"]] <- st2@active.ident

#Save Seurat Object:
save(st2, file = 'IL33_OT3_st2_cells.RData')

#### Preprocessing - Fibroblasts ------------------------------------------------------------------------------------------------------------- ####
#Subset and Re-Cluster Fibroblast and Neighboring Clusters Alone:
Idents(IL33_OT) <- "manual_clusters"
fibroblasts <- subset(IL33_OT, idents = c("Fibroblast", "Mesothelial", "EMT"))

#Find Variable Genes:
fibroblasts <- FindVariableFeatures(fibroblasts, selection.method = "vst", nfeatures = 2000)

#Scale Data:
all.genes <- rownames(fibroblasts)
fibroblasts <- ScaleData(fibroblasts, verbose = T, features = all.genes)

#Run PCA and Determine Dimensions for 90% Variance:
fibroblasts <- RunPCA(fibroblasts, npcs = 30, verbose = FALSE)
st_dev <- fibroblasts@reductions$pca@stdev
var <- st_dev^2
sum(var[1:21])/ sum(var)

#Find Neighbors and Cluster Cells:
fibroblasts <- FindNeighbors(object = fibroblasts, dims = 1:21)
fibroblasts <- FindClusters(object = fibroblasts, resolution = 5)

#Run UMAP:
fibroblasts <- RunUMAP(fibroblasts, reduction = "pca", dims = 1:21, verbose = F)
DimPlot(fibroblasts, reduction = "umap", label = T)

#Identify Cell Populations:
DotPlot(fibroblasts, 
        features = c("Krt19", "Krt8","Krt18","Epcam",
                     "Ero1l", "Vegfa", "Twist1", "Snai1",
                     "Mmp7", "Clu", "Spp1", 
                     "Try4", "Amy2a2", 
                     "Pecam1", "Cdh5",
                     "Msln", "Lrrn4",
                     "Col1a1","Pdpn", "Pdgfra", "Pdgfrb",
                     "Cspg4",
                     "Ptprc",
                     "Cd68","Itgam", "Cd14", "Adgre1", "Mrc1",
                     "Ccr2","Csf1r", "Tlr2",
                     "H2-Eb1","Itgae","Clec9a","Batf3",
                     "S100a8", "Cd33", 
                     "Cd3e", "Cd8a",
                     "Nkg7", "Gzmb",
                     "Cd4",
                     "Gata3", "Il1rl1",
                     "Mitf", "Kit",
                     "Cd79a", "Cd19", "Ms4a1","Ccr10","Prdm1",
                     "Mki67", "Ccna2",
                     "Hbb-bt", "Hbb-bs"), cols = "RdYlBu", dot.scale = 8) + RotatedAxis()

#Remove Mesothelial Cells, EMT-Like cells, and Contaminating Immune Cells:
fibroblasts <- subset(fibroblasts, idents = c(0,2,3,5,8,22,24,25,26), invert = T) 

#Find Variable Genes:
fibroblasts <- FindVariableFeatures(fibroblasts, selection.method = "vst", nfeatures = 2000)

#Scale Data:
all.genes <- rownames(fibroblasts)
fibroblasts <- ScaleData(fibroblasts, verbose = T, features = all.genes)

#Run PCA and Determine Dimensions for 90% Variance:
fibroblasts <- RunPCA(fibroblasts, npcs = 30, verbose = FALSE)
st_dev <- fibroblasts@reductions$pca@stdev
var <- st_dev^2
sum(var[1:22])/ sum(var)

#Find Neighbors and Cluster Cells:
fibroblasts <- FindNeighbors(object = fibroblasts, dims = 1:22)
fibroblasts <- FindClusters(object = fibroblasts, resolution = 4.8)

#Run UMAP:
fibroblasts <- RunUMAP(fibroblasts, reduction = "pca", dims = 1:22, verbose = F)
DimPlot(fibroblasts, reduction = "umap", label = T)

#Identify Cell Populations:
DotPlot(fibroblasts, 
        features = c("Krt19", "Krt8","Krt18","Epcam",
                     "Ero1l", "Vegfa", "Twist1", "Snai1",
                     "Mmp7", "Clu", "Spp1", 
                     "Try4", "Amy2a2", 
                     "Pecam1", "Cdh5",
                     "Msln", "Lrrn4",
                     "Col1a1","Pdpn", "Pdgfra", "Pdgfrb",
                     "Cspg4",
                     "Ptprc",
                     "Cd68","Itgam", "Cd14", "Adgre1", "Mrc1",
                     "Ccr2","Csf1r", "Tlr2",
                     "H2-Eb1","Itgae","Clec9a","Batf3",
                     "S100a8", "Cd33", 
                     "Cd3e", "Cd8a",
                     "Nkg7", "Gzmb",
                     "Cd4",
                     "Gata3", "Il1rl1",
                     "Mitf", "Kit",
                     "Cd79a", "Cd19", "Ms4a1","Ccr10","Prdm1",
                     "Mki67", "Ccna2",
                     "Hbb-bt", "Hbb-bs"), cols = "RdYlBu", dot.scale = 8) + RotatedAxis()

#Remove Proliferating Cells and Non-Fibroblasts:
fibroblasts <- subset(fibroblasts, idents = c(11,20), invert = T) 

#Find Variable Genes:
fibroblasts <- FindVariableFeatures(fibroblasts, selection.method = "vst", nfeatures = 2000)

#Scale Data:
all.genes <- rownames(fibroblasts)
fibroblasts <- ScaleData(fibroblasts, verbose = T, features = all.genes)

#Run PCA and Determine Dimensions for 90% Variance:
fibroblasts <- RunPCA(fibroblasts, npcs = 30, verbose = FALSE)
st_dev <- fibroblasts@reductions$pca@stdev
var <- st_dev^2
sum(var[1:22])/ sum(var)

#Find Neighbors and Cluster Cells:
fibroblasts <- FindNeighbors(object = fibroblasts, dims = 1:22)
fibroblasts <- FindClusters(object = fibroblasts, resolution = 5)

#Run UMAP:
fibroblasts <- RunUMAP(fibroblasts, reduction = "pca", dims = 1:22, verbose = F)
DimPlot(fibroblasts, reduction = "umap", label = T)
DimPlot(fibroblasts, reduction = "umap", label = T, split.by = "Group")

#Identify CAF Populations:
DotPlot(fibroblasts, features = c("Acta2", "Tagln","Clec3b", "Col14a1", "Has1", "Il6","H2-Ab1", "Saa3", "Slpi", "Cd74", "Mki67", "Ccna2"), cols = "RdBu", dot.scale = 8) + RotatedAxis()
VlnPlot(fibroblasts, features = c("Acta2", "Tagln", "Clec3b", "Col14a1", "Has1", "H2-Ab1", "Cd74", "Slpi", "Saa3"), ncol = 3)
VlnPlot(fibroblasts, features = c("Clec3b", "Ly6c1", "Col14a1",
                                  "Acta2", "Col8a1", "Ccn2"), ncol = 3)
VlnPlot(fibroblasts, features = c("Acta2", "Tagln", "Clec3b", "Col14a1"), ncol = 2)

#Label Fibroblast Clusters:
Idents(fibroblasts) <- "seurat_clusters"
fibroblasts <- RenameIdents(fibroblasts,
                            "0" = "apCAF", 
                            "1" = "myCAF", 
                            "2" = "iCAF-hi, myCAF-lo", 
                            "3" = "myCAF-hi, apCAF-lo", 
                            "4" = "iCAF", 
                            "5" = "iCAF-hi, myCAF-lo", 
                            "6" = "myCAF-hi, apCAF-lo", 
                            "7" = "iCAF-hi, myCAF-lo", 
                            "8" = "myCAF-hi, iCAF-lo",
                            "9" = "iCAF-hi, myCAF-lo", 
                            "10" = "myCAF-hi, iCAF-lo",
                            "11" = "myCAF-hi, apCAF-lo",
                            "12" = "iCAF", 
                            "13" = "iCAF", 
                            "14" = "myCAF-hi, iCAF-lo", 
                            "15" = "iCAF", 
                            "16" = "myCAF-hi, apCAF-lo", 
                            "17" = "myCAF-hi, iCAF-lo",
                            "18" = "myCAF-hi, iCAF-lo", 
                            "19" = "iCAF-hi, myCAF-lo", 
                            "20" = "myCAF-hi, iCAF-lo",
                            "21" = "iCAF-hi, myCAF-lo",
                            "22" = "iCAF-hi, myCAF-lo",
                            "23" = "myCAF-hi, iCAF-lo")

new_order <- c("iCAF",
               "iCAF-hi, myCAF-lo",
               "myCAF-hi, iCAF-lo",
               "myCAF",
               "myCAF-hi, apCAF-lo",
               "apCAF")
fibroblasts@active.ident <- factor(fibroblasts@active.ident, levels = new_order)
fibroblasts[["fibro_clusters"]] <- fibroblasts@active.ident

#Save Seurat Object:
save(fibroblasts, file = 'IL33_OT3_fibroblasts.RData')

#### Preprocessing - Macrophages, MonoMacs, and Granulocytes --------------------------------------------------------------------------------- ####
#Subset and Re-Cluster Myeloid Cells of Interest:
Idents(IL33_OT) <- "paper_clusters"
myeloid <- subset(IL33_OT, idents = c("Macrophage", "MonoMac", "Granulocyte"))

#Find Variable Genes:
myeloid <- FindVariableFeatures(myeloid, selection.method = "vst", nfeatures = 2000)

#Scale Data:
all.genes <- rownames(myeloid)
myeloid <- ScaleData(myeloid, verbose = T, features = all.genes)

#Run PCA and Determine Dimensions for 90% Variance:
myeloid <- RunPCA(myeloid, npcs = 30, verbose = FALSE)
st_dev <- myeloid@reductions$pca@stdev
var <- st_dev^2
sum(var[1:20])/ sum(var)

#Find Neighbors and Cluster Cells:
myeloid <- FindNeighbors(object = myeloid, dims = 1:20)
myeloid <- FindClusters(object = myeloid, resolution = 5)

#Run UMAP:
myeloid <- RunUMAP(myeloid, reduction = "pca", dims = 1:20, verbose = F)
DimPlot(myeloid, reduction = "umap", label = T)

#Identify Cell Populations:
DotPlot(myeloid, 
        features = c("Krt19", "Krt8","Krt18","Epcam",
                     "Ero1l", "Vegfa", "Twist1", "Snai1",
                     "Mmp7", "Clu", "Spp1", 
                     "Try4", "Amy2a2", 
                     "Pecam1", "Cdh5",
                     "Msln", "Lrrn4",
                     "Col1a1","Pdpn", "Pdgfra", "Pdgfrb",
                     "Cspg4",
                     "Ptprc",
                     "Cd68","Itgam", "Cd14", "Adgre1", "Mrc1",
                     "Ccr2","Csf1r", "Tlr2",
                     "H2-Eb1","Itgae","Clec9a","Batf3",
                     "S100a8", "Cd33", 
                     "Cd3e", "Cd8a",
                     "Nkg7", "Gzmb",
                     "Cd4",
                     "Gata3", "Il1rl1",
                     "Mitf", "Kit",
                     "Cd79a", "Cd19", "Ms4a1","Ccr10","Prdm1",
                     "Mki67", "Ccna2",
                     "Hbb-bt", "Hbb-bs"), cols = "RdYlBu", dot.scale = 8) + RotatedAxis()

#Remove Contaminating T Cells:
myeloid <- subset(myeloid, idents = 33, invert = T) 

#Find Variable Genes:
myeloid <- FindVariableFeatures(myeloid, selection.method = "vst", nfeatures = 2000)

#Scale Data:
all.genes <- rownames(myeloid)
myeloid <- ScaleData(myeloid, verbose = T, features = all.genes)

#Run PCA and Determine Dimensions for 90% Variance:
myeloid <- RunPCA(myeloid, npcs = 30, verbose = FALSE)
st_dev <- myeloid@reductions$pca@stdev
var <- st_dev^2
sum(var[1:20])/ sum(var)

#Find Neighbors and Cluster Cells:
myeloid <- FindNeighbors(object = myeloid, dims = 1:20)
myeloid <- FindClusters(object = myeloid, resolution = 4.8)

#Run UMAP:
myeloid <- RunUMAP(myeloid, reduction = "pca", dims = 1:20, verbose = F)
DimPlot(myeloid, reduction = "umap", label = T)

#Identify Cell Populations:
DotPlot(myeloid, 
        features = c("Krt19", "Krt8","Krt18","Epcam",
                     "Ero1l", "Vegfa", "Twist1", "Snai1",
                     "Mmp7", "Clu", "Spp1", 
                     "Try4", "Amy2a2", 
                     "Pecam1", "Cdh5",
                     "Msln", "Lrrn4",
                     "Col1a1","Pdpn", "Pdgfra", "Pdgfrb",
                     "Cspg4",
                     "Ptprc",
                     "Cd68","Itgam", "Cd14", "Adgre1", "Mrc1",
                     "Ccr2","Csf1r", "Tlr2",
                     "H2-Eb1","Itgae","Clec9a","Batf3",
                     "S100a8", "Cd33", 
                     "Cd3e", "Cd8a",
                     "Nkg7", "Gzmb",
                     "Cd4",
                     "Gata3", "Il1rl1",
                     "Mitf", "Kit",
                     "Cd79a", "Cd19", "Ms4a1","Ccr10","Prdm1",
                     "Mki67", "Ccna2",
                     "Hbb-bt", "Hbb-bs"), cols = "RdYlBu", dot.scale = 8) + RotatedAxis()

#Remove Proliferating Cells:
myeloid <- subset(myeloid, idents = c(21,32), invert = T)

#Find Variable Genes:
myeloid <- FindVariableFeatures(myeloid, selection.method = "vst", nfeatures = 2000)

#Scale Data:
all.genes <- rownames(myeloid)
myeloid <- ScaleData(myeloid, verbose = T, features = all.genes)

#Run PCA and Determine Dimensions for 90% Variance:
myeloid <- RunPCA(myeloid, npcs = 30, verbose = FALSE)
st_dev <- myeloid@reductions$pca@stdev
var <- st_dev^2
sum(var[1:19])/ sum(var)

#Find Neighbors and Cluster Cells:
myeloid <- FindNeighbors(object = myeloid, dims = 1:19)
myeloid <- FindClusters(object = myeloid, resolution = 5)

#Run UMAP:
myeloid <- RunUMAP(myeloid, reduction = "pca", dims = 1:19, verbose = F)
DimPlot(myeloid, reduction = "umap", label = T)

#Remove Contaminating Epithelial Cells:
cluster_34 <- FindMarkers(myeloid, ident.1 = 34, only.pos = T)
myeloid <- subset(myeloid, idents = 34, invert = T)

#Find Variable Genes:
myeloid <- FindVariableFeatures(myeloid, selection.method = "vst", nfeatures = 2000)

#Scale Data:
all.genes <- rownames(myeloid)
myeloid <- ScaleData(myeloid, verbose = T, features = all.genes)

#Run PCA and Determine Dimensions for 90% Variance:
myeloid <- RunPCA(myeloid, npcs = 30, verbose = FALSE)
st_dev <- myeloid@reductions$pca@stdev
var <- st_dev^2
sum(var[1:20])/ sum(var)

#Find Neighbors and Cluster Cells:
myeloid <- FindNeighbors(object = myeloid, dims = 1:20)
myeloid <- FindClusters(object = myeloid, resolution = 4.8)

#Run UMAP:
myeloid <- RunUMAP(myeloid, reduction = "pca", dims = 1:20, verbose = F)
DimPlot(myeloid, reduction = "umap", label = T)

#Identify Cell Populations:
#macrophage markers from https://www.cell.com/trends/immunology/fulltext/S1471-4906(22)00094-1
DotPlot(myeloid, features = c("Timd4", "Folr2", "Lyve1"), cols = "RdBu", dot.scale = 8) + RotatedAxis() #tissue resident mac_and_mono
DotPlot(myeloid, features = c("Apoe", "Acp5", "Fabp5", "C1qa", "C1qb", "C1qc", "Trem2"), cols = "RdBu", dot.scale = 8) + RotatedAxis() #lipid-associated TAM
DotPlot(myeloid, features = c("Ccl2", "Ccl7", "Ccl8","Cd274","Ifit1","Ifit2","Ifit3","Ifitm1", "Ifitm3", "Il7r", "Nos2","Rsad2","Stat1", "Tnfsf10","Isg15", "Cxcl9", "Cxcl10"), cols = "RdBu", dot.scale = 8) + RotatedAxis() #IFN-TAM
DotPlot(myeloid, features = c("Arg1", "C1qa","Ccl2","Cd63","Clec4d","Il7r","Spp1","Itga4","Vegfa","Mrc1", "Cd274", "Cx3cr1", "Trem2", "Adgre1", "Ccr2"), cols = "RdBu", dot.scale = 8) + RotatedAxis() #Reg-TAM
DotPlot(myeloid, features = c("Il1b", "Ccl3", "Cxcl1", "Cxcl2", "Cxcl3", "Cxcl5", "Ccl20","Ccl3l1", "Il1rn" ,"G0s2", "Inhba", "Spp1"), cols = "RdBu", dot.scale = 8) + RotatedAxis() #Inflam-TAM
DotPlot(myeloid, features = c("Vegfa", "Spp1", "Arg1", "Adam8", "Bnip3", "Mif", "Slc2a1"), cols = "RdBu", dot.scale = 8) + RotatedAxis() #Angio-TAM
DotPlot(myeloid, features = c("Mki67", "Ccnd1", "Ccna2"), cols = "RdBu", dot.scale = 8) + RotatedAxis() #prolif
DotPlot(myeloid, features = c("Ccl2","Ccl9", "Cd14", "Cd300lf","Cxcl10", "F13a1","Fcn1","Fn1","Ifi205","Ifit2","Ifit3","Tgm2", "Il1r2", "Isg20","Itga4", "Ly6c2", "Mgst1","Vcan", "Plaur", "S100a8","S100a9", "S100a12","Sell", "Tlr2", "Thbs1","Ccr2","Adgre1", "Mrc1"), cols = "RdBu", dot.scale = 8) + RotatedAxis() #Tumor Infiltrating Monocyte
DotPlot(myeloid, features = c("Ace", "Adgre4", "Cd300a", "Cdkn1c", "Ceacam1", "Ear2", "Il17ra", "Itgal", "Lilrb2", "Lrp1", "Spn", "Stk10", "Tnfrsf1b","Trml4", "Ccr2","Adgre1"), cols = "RdBu", dot.scale = 8) + RotatedAxis() #Nonclassical Monoycte

VlnPlot(myeloid, features = c("Ccr2", "Adgre1", "Tnfrsf1b","Trml4"), ncol = 2)

DotPlot(myeloid, features = c("Timd4", "Folr2", "Lyve1", #RTM
                              "Apoe", "Acp5", "Fabp5", "Trem2", #lipid TAM
                              "Arg1", "Cx3cr1", "Mrc1", "Gpnmb", #reg TAM
                              "Cd14", "Ccr2",  #Classical TIM
                              "Ace", "Spn", "Ly6g", #nonclassical mono
                              "Ccna2", "Mki67" #prolif
), cols = "RdYlBu", dot.scale = 8) + RotatedAxis()

FeaturePlot(myeloid, features = "Adgre1", cols = c("gainsboro", "firebrick1"))

#Annotate Myeloid Clusters:
Idents(myeloid) <- "seurat_clusters"
myeloid <- RenameIdents(myeloid,
                        "0" = "Granulocyte", 
                        "1" = "Macrophage", 
                        "2" = "Macrophage", 
                        "3" = "MonoMac", 
                        "4" = "MonoMac", 
                        "5" = "MonoMac", 
                        "6" = "Macrophage", 
                        "7" = "MonoMac", 
                        "8" = "Macrophage",
                        "9" = "Granulocyte", 
                        "10" = "Granulocyte",
                        "11" = "Macrophage",
                        "12" = "Macrophage", 
                        "13" = "Granulocyte", 
                        "14" = "Macrophage", 
                        "15" = "MonoMac", 
                        "16" = "MonoMac", 
                        "17" = "MonoMac",
                        "18" = "Granulocyte", 
                        "19" = "Macrophage", 
                        "20" = "MonoMac",
                        "21" = "Macrophage",
                        "22" = "MonoMac",
                        "23" = "Macrophage",
                        "24" = "Macrophage", 
                        "25" = "MonoMac", 
                        "26" = "MonoMac", 
                        "27" = "Granulocyte",
                        "28" = "MonoMac", 
                        "29" = "MonoMac", 
                        "30" = "Macrophage",
                        "31" = "Macrophage",
                        "32" = "Granulocyte",
                        "33" = "Granulocyte",
                        "34" = "MonoMac", 
                        "35" = "Macrophage",
                        "36" = "Granulocyte", 
                        "37" = "Macrophage")

new_order <- c("Granulocyte",
               "Macrophage",
               "MonoMac")
myeloid@active.ident <- factor(myeloid@active.ident, levels = new_order)
myeloid[["myeloid_clusters"]] <- myeloid@active.ident

#Save Seurat Object:
save(myeloid, file = 'IL33_OT3_mac_mono_granulo.RData')

#### Preprocessing - Tumor Cells ------------------------------------------------------------------------------------------------------------- ####
#Subset and Re-Cluster Tumor Cells:
Idents(IL33_OT) <- "paper_clusters"
tumor <- subset(IL33_OT, idents = c("Tumor", "Proliferating"))

#Find Variable Genes:
tumor <- FindVariableFeatures(tumor, selection.method = "vst", nfeatures = 2000)

#Scale Data:
all.genes <- rownames(tumor)
tumor <- ScaleData(tumor, verbose = T, features = all.genes)

#Run PCA and Determine Dimensions for 90% Variance:
tumor <- RunPCA(tumor, npcs = 30, verbose = FALSE)
st_dev <- tumor@reductions$pca@stdev
var <- st_dev^2
sum(var[1:19])/ sum(var)

#Find Neighbors and Cluster Cells:
tumor <- FindNeighbors(object = tumor, dims = 1:19)
tumor <- FindClusters(object = tumor, resolution = 5)

#Run UMAP:
tumor <- RunUMAP(tumor, reduction = "pca", dims = 1:19, verbose = F)
DimPlot(tumor, reduction = "umap", label = T)

#Identify Cell Populations:
DotPlot(tumor, features = c("Krt19", "Epcam","Msln","Try4", "Col1a2", "Acta2", "Clec3b", "Cspg4","Saa3","Ptprc", "Cd3e", "Trdc","Cd4", "Foxp3", "Ccr7", "Cd8a", "Nkg7","Il2ra","Il1rl1","Arg1", "Ccr2","Cd68","Fcgr3","Adgre1", "Itgam","Cd14","Mrc1", "S100a8", "Cd33","H2-Eb1",
                            "Itgae","Clec9a","Batf3", "Ly6c2", "Ly6g","Cd79a", "Cd19", "Ms4a1","Ccr10","Prdm1","Kit","Ccna2", "Mki67","Pecam1", "Cdh5", "Hbb-bt"), cols = "RdYlBu", dot.scale = 8) + RotatedAxis()

#Remove Non-Tumor Cells:
tumor <- subset(tumor, idents = c(6, 32, 23, 29, 26, 31, 28, 1, 10, 5, 33), invert = T)

#Find Variable Genes:
tumor <- FindVariableFeatures(tumor, selection.method = "vst", nfeatures = 2000)

#Scale Data:
all.genes <- rownames(tumor)
tumor <- ScaleData(tumor, verbose = T, features = all.genes)

#Run PCA and Determine Dimensions for 90% Variance:
tumor <- RunPCA(tumor, npcs = 30, verbose = FALSE)
st_dev <- tumor@reductions$pca@stdev
var <- st_dev^2
sum(var[1:21])/ sum(var)

#Find Neighbors and Cluster Cells:
tumor <- FindNeighbors(object = tumor, dims = 1:21)
tumor <- FindClusters(object = tumor, resolution = 5)

#Run UMAP:
tumor <- RunUMAP(tumor, reduction = "pca", dims = 1:21, verbose = F)
DimPlot(tumor, reduction = "umap", label = T)

#Identify Cell Populations:
DotPlot(tumor, features = c("Krt19", "Epcam","Msln","Try4", "Col1a2", "Acta2", "Clec3b", "Cspg4","Saa3","Ptprc", "Cd3e", "Trdc","Cd4", "Foxp3", "Ccr7", "Cd8a", "Nkg7","Il2ra","Il1rl1","Arg1", "Ccr2","Cd68","Fcgr3","Adgre1", "Itgam","Cd14","Mrc1", "S100a8", "Cd33","H2-Eb1",
                            "Itgae","Clec9a","Batf3", "Ly6c2", "Ly6g","Cd79a", "Cd19", "Ms4a1","Ccr10","Prdm1","Kit","Ccna2", "Mki67","Pecam1", "Cdh5", "Hbb-bt"), cols = "RdYlBu", dot.scale = 8) + RotatedAxis()

#Remove Non-Tumor Cells:
tumor <- subset(tumor, idents = c(27,28), invert = T) #includes proliferating cells

#Save Seurat Object:
save(tumor, file = 'IL33_OT3_tumor.RData')

#### Preprocessing - CellChat Interaction Mapping -------------------------------------------------------------------------------------------- ####
Idents(object = IL33_OT) <- "paper_clusters"

#Add Treg Identities from ST2 Object:
treg_ids <- WhichCells(st2, idents = "Treg")
Idents(IL33_OT, cells = treg_ids) <- "Treg"

#Remove Proliferating Cells and Red Blood Cells:
cellchat <- subset(IL33_OT, idents = c("Proliferating", "RBC"), invert = T)

#Create CellChat Cluster Annotation:
new_order <- c("Tumor",
               "EMT",
               "Ductal",
               "Acinar",
               "Endothelial",
               "Mesothelial",
               "Fibroblast",
               "Pericyte",
               "Macrophage",
               "MonoMac",
               "Dendritic Cell",
               "Granulocyte",
               "CD8+ T Cell",
               "NK Cell",
               "CD4+ T Cell",
               "Treg",
               "ILC2",
               "Mast Cell",
               "B Cell",
               "Plasma Cell")
cellchat@active.ident <- factor(cellchat@active.ident, levels = new_order)
cellchat[["cc_clusters"]] <- cellchat@active.ident

#Create CellChat Objects:
Idents(object = cellchat) <- "Group"
CC_CTL <- subset(x = cellchat, idents ="Pdgfra-CreERT2")
CC_Il33 <- subset(x = cellchat, idents ="Pdgfra-CreERT2;IL33 f/f")

cellchat_CTRL <- createCellChat(object = CC_CTL, group.by = "cc_clusters", assay = "RNA")
cellchat_IL33 <- createCellChat(object = CC_Il33, group.by = "cc_clusters", assay = "RNA")

#Import Ligand-Receptor Interaction Database:
CellChatDB <- CellChatDB.mouse
cellchat_CTRL@DB <- CellChatDB
cellchat_IL33@DB <- CellChatDB

#Preprocess Expression Data:
cellchat_CTRL <- subsetData(cellchat_CTRL) 
future::plan("multisession", workers = 4) 

cellchat_IL33 <- subsetData(cellchat_IL33) 
future::plan("multisession", workers = 4) 

cellchat_CTRL <- identifyOverExpressedGenes(cellchat_CTRL)
cellchat_CTRL <- identifyOverExpressedInteractions(cellchat_CTRL)

cellchat_IL33 <- identifyOverExpressedGenes(cellchat_IL33)
cellchat_IL33 <- identifyOverExpressedInteractions(cellchat_IL33)

#Compute Communication Probability and Infer Cellular Communication Network:
cellchat_CTRL <- computeCommunProb(cellchat_CTRL, population.size = T)
cellchat_CTRL <- filterCommunication(cellchat_CTRL, min.cells = 5)

cellchat_IL33 <- computeCommunProb(cellchat_IL33, population.size = T)
cellchat_IL33 <- filterCommunication(cellchat_IL33, min.cells = 5)

#Compute Communication Probability at Pathway Level:
cellchat_CTRL <- computeCommunProbPathway(cellchat_CTRL)
cellchat_IL33 <- computeCommunProbPathway(cellchat_IL33)

#Compute Total Number of Interactions:
cellchat_CTRL <- aggregateNet(cellchat_CTRL)
cellchat_IL33 <- aggregateNet(cellchat_IL33)

cellchat_CTRL <- netAnalysis_computeCentrality(cellchat_CTRL)
cellchat_IL33 <- netAnalysis_computeCentrality(cellchat_IL33)

#Merge the Two CellChat Objects
object.list <- list(CTRL = cellchat_CTRL, IL33 = cellchat_IL33)
cellchat_MERGE <- mergeCellChat(object.list, add.names = names(object.list))
saveRDS(cellchat_MERGE, file = "cellchat_comparison_merge_OT3.rds")
saveRDS(object.list, file = "cellchat_object.list_OT3.rds")

#Perform DE Ligand Analysis Between Groups
pos.dataset = "IL33"
features.name = pos.dataset
cellchat_MERGE <- identifyOverExpressedGenes(cellchat_MERGE, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.01, thresh.fc = 0.01, thresh.p = 0.05)

#Map the Differential Expression Analysis Results
net <- netMappingDEG(cellchat_MERGE, features.name = features.name)

#Extract the Ligand-Receptor Pairs with Upregulated Ligands in IL33
net.up <- subsetCommunication(cellchat_MERGE, net = net, datasets = "IL33",ligand.logFC = 0.25, receptor.logFC = NULL)
gene.up <- extractGeneSubsetFromPair(net.up, cellchat_MERGE)

#Extract the Ligand-Receptor Pairs with Upregulated Ligands in CTRL
net.down <- subsetCommunication(cellchat_MERGE, net = net, datasets = "CTRL",ligand.logFC = -0.25, receptor.logFC = NULL)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat_MERGE)

#### Figures---------------------------------------------------------------------------------------------------------------------------------- ####
global_colors <- c("Tumor" = "#ffd92f",
                   "EMT" = "#bf9b30",
                   "Acinar" = "#a66999",
                   "Ductal" = "#D82C20",
                   "Endothelial" = "#6b91c9",
                   "Mesothelial" = "#006666",
                   "Fibroblast" = "#66b2b2",
                   "Pericyte" = "#66ccff",
                   "Macrophage" = "#fdb462",
                   "Dendritic Cell" = "#f6e0b5",
                   "MonoMac" = "#cc6633",
                   "Granulocyte" = "#decbe4",
                   "CD4+ T Cell" = "#b3de69",
                   "ILC2" = "#4daf4a",
                   "CD8+ T Cell" = "#fccde5",
                   "NK Cell" = "#ef008c",
                   "Mast Cell" = "#0000ff",
                   "B Cell" = "#66545e",
                   "Plasma Cell" = "#a39193",
                   "Proliferating" = "#dfe3ee",
                   "RBC" = "#dec2cb")

plot_colors <- c("Tumor" = "#ffd92f",
                 "EMT" = "#bf9b30",
                 "Acinar" = "#a66999",
                 "Ductal" = "#D82C20",
                 "Endothelial" = "#6b91c9",
                 "Mesothelial" = "#006666",
                 "Fibroblast" = "#66b2b2",
                 "Pericyte" = "#66ccff",
                 "Macrophage" = "#fdb462",
                 "Dendritic Cell" = "#f6e0b5",
                 "MonoMac" = "#cc6633",
                 "Granulocyte" = "#decbe4",
                 "CD4+ T Cell" = "#b3de69",
                 "Treg" = "maroon",
                 "ILC2" = "#4daf4a",
                 "CD8+ T Cell" = "#fccde5",
                 "NK Cell" = "#ef008c",
                 "Mast Cell" = "#0000ff",
                 "B Cell" = "#66545e",
                 "Plasma Cell" = "#a39193")

#### Figure 2H, Split Global UMAPs #### 
Idents(object = IL33_OT) <- 'paper_clusters'
DimPlot(IL33_OT, label = F, cols = global_colors, pt.size = 0.2,order = rev(c("Tumor",
                                                                              "Ductal",
                                                                              "Endothelial",
                                                                              "Fibroblast",
                                                                              "Mesothelial",
                                                                              "EMT",
                                                                              "Pericyte",
                                                                              "MonoMac",
                                                                              "Macrophage",
                                                                              "Acinar",
                                                                              "Dendritic Cell",
                                                                              "Mast Cell",
                                                                              "Granulocyte",
                                                                              "CD4+ T Cell",
                                                                              "CD8+ T Cell",
                                                                              "NK Cell",
                                                                              "ILC2",
                                                                              "B Cell",
                                                                              "Plasma Cell",
                                                                              "Proliferating",
                                                                              "RBC")), split.by = "Group")
#### Figure 2I, Tumor Cell Pathway Analysis ####
Idents(IL33_OT) <- "paper_clusters"
tumor <- subset(IL33_OT, idents = "Tumor")
tumor_DE <- FindMarkers(tumor, ident.1 = "Pdgfra-CreERT2", ident.2 = "Pdgfra-CreERT2;IL33 f/f", group.by = "Group", logfc.threshold = 0, min.pct = 0)

h_gene_sets = msigdbr(species = "mouse", category = "H")
msigdbr_list = split(x = h_gene_sets$gene_symbol, f = h_gene_sets$gs_name)

ranks <- tumor_DE$avg_log2FC
names(ranks) <- rownames(tumor_DE)

fgseaRes <- fgsea(pathways = msigdbr_list, 
                  stats = ranks)

result <- apply(fgseaRes,2,as.character)
write.csv(result, file = "OT3 Tumor Cell CreERT2 VS CreERT2;Il33 fGSEA Result.csv")

fgseaRes_sig <- subset(fgseaRes, fgseaRes$padj <= 0.05)

waterfall_plot_sig <- function (fgseaRes, graph_title) {
  (fgseaRes %>% 
     mutate(short_name = str_split_fixed(pathway, "_",2)[,2])%>%
     ggplot( aes(reorder(short_name,NES), NES)) +
     geom_bar(stat= "identity", aes(fill = padj))+
     coord_flip()+
     labs(x = "Hallmark Gene Set", y = "Normalized Enrichment Score", title = graph_title)+
     theme(axis.text.y = element_text(size = 7), 
           plot.title = element_text(hjust = 1)))+scale_fill_viridis(direction = 1, option = "A")
}
waterfall_plot_sig(fgseaRes_sig, "Pathways Enriched in WT vs. KO - Tumor Cells")

#### Figure 2J, Il1rl1 Violin Plots ####
VlnPlot(st2, features = "Il1rl1", split.by = "Group", cols = c("#dfeefb", "#5fa8e9"))

#### Figure S2D, Global Cluster Marker Dotplot ####
Idents(IL33_OT) <- "paper_clusters"
DotPlot(IL33_OT, 
        features = c("Krt19", "Krt8","Krt18","Epcam",
                     "Ero1l", "Vegfa", "Twist1", "Snai1",
                     "Mmp7", "Clu", "Spp1", 
                     "Try4", "Amy2a2", 
                     "Pecam1", "Cdh5",
                     "Msln", "Lrrn4",
                     "Col1a1","Pdpn", "Pdgfra", "Pdgfrb",
                     "Cspg4",
                     "Ptprc",
                     "Cd68", "Msr1", "C1qc", "Mrc1", "Adgre1", "Sirpa", "Fcgr4", "Csf1r", "Lgals3",
                     "Itgam", "Cd14","Ccr2", 
                     "H2-Eb1","Itgae","Clec9a","Batf3",
                     "S100a8", "Cd33", 
                     "Cd3e", "Cd8a",
                     "Nkg7", "Gzmb",
                     "Cd4",
                     "Gata3", "Il1rl1",
                     "Mitf", "Kit",
                     "Cd79a", "Cd19", "Ms4a1","Ccr10","Prdm1",
                     "Mki67", "Ccna2",
                     "Hbb-bt", "Hbb-bs"), cols = "RdYlBu", dot.scale = 8) + RotatedAxis()

#### Figure S2E, Cell Abundance Bar Graph ####
Idents(object = IL33_OT) <- 'Group'
samples.list <- unique(IL33_OT$Group)
clusters <- lapply(samples.list, function(x){
  subset <- subset(IL33_OT, subset = Group == x)
  dist <- data.frame(table(subset$paper_clusters))
  
  return(dist)
})

names(clusters) <- samples.list

clusters_percent <- lapply(clusters, FUN = function(x){
  summ <- sum(x$Freq)
  x$Freq <- (x$Freq/summ)
  return(x)
})

clusters_dist <- reshape2::melt(clusters, id.var = "Var1")
colnames(clusters_dist) <- c("paper_clusters","variable","value","Group")
clusters_percent_dist <- reshape2::melt(clusters_percent, id.var = "Var1")
colnames(clusters_percent_dist) <- c("paper_clusters","variable","value","Group")

ggplot(clusters_percent_dist, aes(fill=paper_clusters, y = value, x = Group)) + 
  scale_x_discrete(limits = c("Pdgfra-CreERT2", "Pdgfra-CreERT2;IL33 f/f")) + 
  geom_bar(position = "stack", stat = "identity") + 
  theme_bw()+
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, size = 8, face = "bold", hjust = 1))+
  ylab("count") + xlab("orig.ident") + ggtitle("Relative Cell Types Abundance") +scale_fill_manual(values = global_colors)


#### Figure S2F, Global Il33, Pdgfra, Pdgfrb Violin Plots ####
Stacked_VlnPlot(seurat_object = IL33_OT, features = c("Il33", "Pdgfra", "Pdgfrb"), x_lab_rotate = TRUE,
                colors_use =c("#dfeefb", "#5fa8e9") , split.by = "Group")

fibroblasts <- subset(IL33_OT, idents = "Fibroblast")
mesothelial <- subset(IL33_OT, idents = "Mesothelial")
tumor <- subset(IL33_OT, idents = "Tumor")
peri <- subset(IL33_OT, idents = "Pericyte")
emt <- subset(IL33_OT, idents = "EMT")

fibro_de <- FindMarkers(fibroblasts, ident.1 = "Pdgfra-CreERT2", ident.2 = "Pdgfra-CreERT2;IL33 f/f", group.by = "Group" )
meso_de <- FindMarkers(mesothelial, ident.1 = "Pdgfra-CreERT2", ident.2 = "Pdgfra-CreERT2;IL33 f/f", group.by = "Group" )
tumor_de <- FindMarkers(tumor, ident.1 = "Pdgfra-CreERT2", ident.2 = "Pdgfra-CreERT2;IL33 f/f", group.by = "Group" )
peri_de <- FindMarkers(peri, ident.1 = "Pdgfra-CreERT2", ident.2 = "Pdgfra-CreERT2;IL33 f/f", group.by = "Group" )
emt_de <- FindMarkers(emt, ident.1 = "Pdgfra-CreERT2", ident.2 = "Pdgfra-CreERT2;IL33 f/f", group.by = "Group" )

fibro_de["genes"]<- row.names(fibro_de)
#Il33 4.205107e-61
#pdgfra ns
#pdgfrb ns

meso_de["genes"]<- row.names(meso_de)
#Il33 1.588767e-02
#pdgfra ns
#pdgfrb ns

tumor_de["genes"]<- row.names(tumor_de)
#Il33 ns
#pdgfra ns
#pdgfrb ns

peri_de["genes"]<- row.names(peri_de)
#Il33 ns
#pdgfra ns
#pdgfrb ns

emt_de["genes"]<- row.names(emt_de)
#Il33 0.03019773
#pdgfra ns
#pdgfrb ns

#### Figure S2G, Xist Violin Plot ####
xist_cells <- subset(IL33_OT, idents = c("Tumor", "EMT", "Ductal", "Fibroblast", "Mesothelial", "Pericyte"))
Stacked_VlnPlot(xist_cells, features = "Xist", colors_use = plot_colors)

#### Figure S2H, Classical and Basal Tumor Cell Scoring ####
basal_sig <- c("S100a2", "Ly6d", "Sprr1b", 
               "Lemd1", "Krt15", "Ctsl", "Dhrs9", "Areg", "Cst6", 
               "Fam83a", "Scel", "Fgfbp1", "Krt7",
               "Krt17", "Gpr87", "Tns4", "Slc2a1", "Anxa8")

classical_sig <- c("Oit1", "Prr15l", "Ctse", "Lyz1", "Lyz2",
                   "Tff2", "Tff1", "Anxa10", "Lgals4", "Pla2g10", "Ceacam1","Ceacam2", "Psg16", "Vsig2",
                   "Tspan8", "St6galnac1", "Agr2", "Tff3", "Myo1a", 
                   "Clrn3", "Krt20", "Cdh17", "Spink4", "Reg4",
                   "Cyp3a13","Cyp3a57")

Idents(tumor) <- "Group"
levels(tumor@active.ident)

#Add Module Scores to New Object:
tumor_2 <- AddModuleScore(tumor,
                          features = list(classical_sig),
                          name="Classical_Score")
tumor_2 <- AddModuleScore(tumor_2,
                          features = list(basal_sig),
                          name="Basal_Score")

#Visualize Module Scores:
pal <- brewer.pal(9, "PuRd")
FeaturePlot_scCustom(seurat_object = tumor_2, features = c("Classical_Score1", "Basal_Score1"), na_cutoff = NULL, order = T, split.by = "Group", colors_use = pal)

#Calculate p values:
score_data <- FetchData(tumor_2, vars = c("Classical_Score1","Basal_Score1", "Group"))

p_value <- pairwise.wilcox.test(score_data$Classical_Score1, score_data$Group, p.adjust.method="bonferroni")
p_value$p.value
#2.501504e-09

p_value2 <- pairwise.wilcox.test(score_data$Basal_Score1, score_data$Group, p.adjust.method="bonferroni")
p_value2$p.value
#0.4924549

#### Figure S3D, ST2 Cluster Marker Dotplot ####
DotPlot(st2, features = c("Ptprc", "Cd3e", "Cd4", 
                          "Tcf7", "Ccr7", "Sell",
                          "Pdcd1", "Tox", "Cd160",
                          "Tbx21", "Cxcr3", "Ifng",
                          "Rorc", "Ccr6", "Il17a", 
                          "Foxp3", "Il2ra",
                          "Gata3", "Kit", "Mitf"), cols = "RdBu", dot.scale = 9) + RotatedAxis()

#### Figure 3A, ST2 Activation Dotplot ####
st2_only <- subset(st2, idents = c("Treg", "ILC2", "Mast Cell")) 
new_order <- c("ILC2",
               "Mast Cell",
               "Treg")
st2_only@active.ident <- factor(st2_only@active.ident, levels = new_order)
st2_only[["st2_clusters"]] <- st2_only@active.ident
DotPlot(st2_only, features = c("Il4", "Il5", "Il13", "Lif", "Csf2", "Ccl3","Il6","Areg", "Klrg1"), cols = "RdYlBu", dot.scale = 10, split.by = 'Group') + RotatedAxis()

#### Figure 3B, Fibroblast Pathway Analysis ####
fb_DE <- FindMarkers(fibroblasts, ident.1 = "Pdgfra-CreERT2", ident.2 = "Pdgfra-CreERT2;IL33 f/f", group.by = "Group", logfc.threshold = 0, min.pct = 0)

ranks <- fb_DE$avg_log2FC
names(ranks) <- rownames(fb_DE)

h_gene_sets = msigdbr(species = "mouse", category = "H")
msigdbr_list = split(x = h_gene_sets$gene_symbol, f = h_gene_sets$gs_name)

fgseaRes <- fgsea(pathways = msigdbr_list, 
                  stats = ranks)

result <- apply(fgseaRes,2,as.character)
write.csv(result, file = "OT3 Fibroblast CreERT2 VS CreERT2;Il33 fGSEA Result HALLMARK.csv")

fgseaRes_sig <- subset(fgseaRes, fgseaRes$padj <= 0.05)

waterfall_plot_sig <- function (fgseaRes, graph_title) {
  (fgseaRes %>% 
     mutate(short_name = str_split_fixed(pathway, "_",2)[,2])%>%
     ggplot( aes(reorder(short_name,NES), NES)) +
     geom_bar(stat= "identity", aes(fill = padj))+
     coord_flip()+
     labs(x = "Hallmark Gene Set", y = "Normalized Enrichment Score", title = graph_title)+
     theme(axis.text.y = element_text(size = 7), 
           plot.title = element_text(hjust = 1)))+scale_fill_viridis(direction = 1, option = "G")
}
waterfall_plot_sig(fgseaRes_sig, "Pathways Enriched in WT vs. KO - Fibroblasts")

#### Figure 3C, Fibroblast Interaction Chord Diagrams - Receiving ####
strwidth <- function(x) {0.01}
netVisual_chord_gene(object.list[[1]], sources.use = c(1:20), targets.use = 7, slot.name = 'net', net = net.down, scale = F, link.target.prop = T, lab.cex = 0.5,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]), color.use = plot_colors, show.legend = F)
netVisual_chord_gene(object.list[[1]], sources.use = c(1:20), targets.use = 7, slot.name = 'netP', net = net.down, scale = F, link.target.prop = T, lab.cex = 0.5,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]), color.use = plot_colors, show.legend = F)

netVisual_chord_gene(object.list[[2]], sources.use = c(1:20), targets.use = 7, slot.name = 'net', net = net.up, scale = F, link.target.prop = T,lab.cex = 0.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]), color.use = plot_colors, show.legend = F)
netVisual_chord_gene(object.list[[2]], sources.use = c(1:20), targets.use = 7, slot.name = 'netP', net = net.up, scale = F, link.target.prop = T,lab.cex = 0.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]), color.use = plot_colors, show.legend = F)


#### Figure 3D, Fibroblast UMAP ####
DimPlot(fibroblasts, reduction = "umap", label = F,
        cols = c("iCAF" = "#efa14b",
                 "iCAF-hi, myCAF-lo" = "#ffd6a8",
                 "myCAF-hi, iCAF-lo" = "#91d0ff",
                 "myCAF" = '#2986CC',
                 "myCAF-hi, apCAF-lo" = "#9dd4d4",
                 "apCAF" = "#099393"), group.by = 'fibro_clusters')

#### Figure 3E, CAF Marker Violin Plots ####
Idents(fibroblasts) <- "fibro_clusters"
Stacked_VlnPlot(seurat_object = fibroblasts, features = c("Pdpn", "Pdgfra","Pdgfrb",
                                                          "Clec3b", "Ly6c1", "Col14a1"), colors_use = c("iCAF" = "#efa14b",
                                                                                                        "iCAF-hi, myCAF-lo" = "#ffd6a8",
                                                                                                        "myCAF-hi, iCAF-lo" = "#91d0ff",
                                                                                                        "myCAF" = '#2986CC',
                                                                                                        "myCAF-hi, apCAF-lo" = "#9dd4d4",
                                                                                                        "apCAF" = "#099393"), 
                x_lab_rotate = TRUE)

Stacked_VlnPlot(seurat_object = fibroblasts, features = c("Acta2", "Col8a1", "Ccn2",
                                                          "Saa3", "H2-Ab1", "Slpi"), colors_use = c("iCAF" = "#efa14b",
                                                                                                    "iCAF-hi, myCAF-lo" = "#ffd6a8",
                                                                                                    "myCAF-hi, iCAF-lo" = "#91d0ff",
                                                                                                    "myCAF" = '#2986CC',
                                                                                                    "myCAF-hi, apCAF-lo" = "#9dd4d4",
                                                                                                    "apCAF" = "#099393"), 
                x_lab_rotate = TRUE)

#### Figure 3F, Il33 Split Violin Plots ####
Idents(fibroblasts) <- "fibro_clusters"
Stacked_VlnPlot(seurat_object = fibroblasts, features = "Il33", x_lab_rotate = TRUE,
                colors_use =c("#dfeefb", "#5fa8e9") , split.by = "Group")

#### Figure 3G, Fibroblast Cell Abundance Bar Graph ####
Idents(object = fibroblasts) <- 'Group'
samples.list <- unique(fibroblasts$Group)
clusters <- lapply(samples.list, function(x){
  subset <- subset(fibroblasts, subset = Group == x)
  dist <- data.frame(table(subset$fibro_clusters))
  
  return(dist)
})

names(clusters) <- samples.list

clusters_percent <- lapply(clusters, FUN = function(x){
  summ <- sum(x$Freq)
  x$Freq <- (x$Freq/summ)
  return(x)
})

clusters_dist <- reshape2::melt(clusters, id.var = "Var1")
colnames(clusters_dist) <- c("fibro_clusters","variable","value","Group")
clusters_percent_dist <- reshape2::melt(clusters_percent, id.var = "Var1")
colnames(clusters_percent_dist) <- c("fibro_clusters","variable","value","Group")

ggplot(clusters_percent_dist, aes(fill=fibro_clusters, y = value, x = Group)) + 
  scale_x_discrete(limits = c("Pdgfra-CreERT2", "Pdgfra-CreERT2;IL33 f/f")) + 
  geom_bar(position = "stack", stat = "identity") + 
  theme_bw()+
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, size = 8, face = "bold", hjust = 1))+
  ylab("count") + xlab("orig.ident") + ggtitle("Relative Cell Types Abundance") +scale_fill_manual(values = c("iCAF" = "#efa14b",
                                                                                                              "iCAF-hi, myCAF-lo" = "#ffd6a8",
                                                                                                              "myCAF-hi, iCAF-lo" = "#91d0ff",
                                                                                                              "myCAF" = '#2986CC',
                                                                                                              "myCAF-hi, apCAF-lo" = "#9dd4d4",
                                                                                                              "apCAF" = "#099393"))

#### Figure 3H, Fibroblast Interaction Chord Diagrams - Sending ####
strwidth <- function(x) {0.01}

netVisual_chord_gene(object.list[[1]], sources.use = 7, targets.use = c(1:20), slot.name = 'net', net = net.down, scale = F, link.target.prop = T, lab.cex = 0.5,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]), color.use = plot_colors, show.legend = F)
netVisual_chord_gene(object.list[[1]], sources.use = 7, targets.use = c(1:20), slot.name = 'netP', net = net.down, scale = F, link.target.prop = T, lab.cex = 0.5,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]), color.use = plot_colors, show.legend = F)

netVisual_chord_gene(object.list[[2]], sources.use = 7, targets.use = c(1:20), slot.name = 'net', net = net.up, scale = F, link.target.prop = T,lab.cex = 0.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]), color.use = plot_colors, show.legend = F)
netVisual_chord_gene(object.list[[2]], sources.use = 7, targets.use = c(1:20), slot.name = 'netP', net = net.up, scale = F, link.target.prop = T,lab.cex = 0.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]), color.use = plot_colors, show.legend = F)

#### Figure S4A, Egfr and Areg Violin Plots ####
Idents(cellchat) <- "cc_clusters"
Stacked_VlnPlot(seurat_object = cellchat, features = c("Areg", "Egfr"), x_lab_rotate = TRUE,
                colors_use =c("#dfeefb", "#5fa8e9") , split.by = "Group")

#### Figure S4G, CAF Ligand Dotplot ####
DotPlot(fibroblasts, features = c("Col1a1", "Col1a2","C3", "Ptn", "Tnxb", "Postn","Cxcl5", "Cxcl12", "Igf1", "Vegfd", "Clec2d", "Gdf10", "Kng2",
                                  "Thbs1", "Mif", "Spp1", "Tnc", "Ncam1", "Cxcl1", "Cxcl2", "Vegfa", "Cx3cl1", "Tgfb1", "Il1b"), cols = "RdYlBu", dot.scale = 8, split.by = "Group") + RotatedAxis()

#### Figure 4B, Macrophage Activation Dotplot ####
mac_and_mono <- subset(myeloid, idents = c("Macrophage", "MonoMac"))
mac_and_mono <- ScaleData(mac_and_mono, verbose = T, features = row.names(mac_and_mono))
DotPlot(mac_and_mono, features = c("Il1a", "Il1b","Cxcl10", "Cd80", "Cd86","Tnf","Tlr2", "Il6","Cxcl9", "H2-Ab1","H2-Eb1", "Ccl5","Cxcl16","Tlr4",     
                                   "Tgfb1", "Vegfa",  "Ccl24", "Chil3", "Arg1", "Mrc1",  "Cd163"), cols = "RdBu", dot.scale = 10, split.by = "Group") + RotatedAxis()

#### Figure 4C, CD8 T Cell Interaction Chord Diagrams ####
strwidth <- function(x) {0.01}

netVisual_chord_gene(object.list[[1]], sources.use = c(1:20), targets.use = 13, slot.name = 'net', net = net.down, scale = F, link.target.prop = T, lab.cex = 0.5,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]), color.use = plot_colors, show.legend = F)
netVisual_chord_gene(object.list[[1]], sources.use = c(1:20), targets.use = 13, slot.name = 'netP', net = net.down, scale = F, link.target.prop = T, lab.cex = 0.5,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]), color.use = plot_colors, show.legend = F)

netVisual_chord_gene(object.list[[2]], sources.use = c(1:20), targets.use = 13, slot.name = 'net', net = net.up, scale = F, link.target.prop = T,lab.cex = 0.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]), color.use = plot_colors, show.legend = F)
netVisual_chord_gene(object.list[[2]], sources.use = c(1:20), targets.use = 13, slot.name = 'netP', net = net.up, scale = F, link.target.prop = T,lab.cex = 0.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]), color.use = plot_colors, show.legend = F)

#### Figure S5B, Granulocyte DE Volcano Plot ####
granu <- subset(myeloid, idents = "Granulocyte")

granu_de <- FindMarkers(granu, ident.1 = "Pdgfra-CreERT2", ident.2 = "Pdgfra-CreERT2;IL33 f/f", group.by = "Group", logfc.threshold = 0, min.pct = 0 )
granu_de_volcano <- as.data.frame(granu_de$p_val_adj)
colnames(granu_de_volcano) <- "p_val_adj"
granu_de_volcano["avg_log2FC"] <- as.data.frame(granu_de$avg_log2FC)
rownames(granu_de_volcano) <- row.names(granu_de)

EnhancedVolcano(granu_de_volcano,
                lab = rownames(granu_de_volcano),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                FCcutoff = 0.5,
                title = 'Granulocyte WT vs KO',
                pCutoff = 0.05,
                ylab = "-log10(adj. p val.)",
                xlim = c(-1.5, 2), 
                ylim = c(0,25),
                labSize = 3.5,
                col = c( "#0065a2","#E7C582", "#C05780",  "#00b0ba"),
                colAlpha = 0.8)

#### Figure S5C, Granulocyte Pathway Analysis ####
h_gene_sets = msigdbr(species = "mouse", category = "H")
msigdbr_list = split(x = h_gene_sets$gene_symbol, f = h_gene_sets$gs_name)

ranks <- granu_de$avg_log2FC
names(ranks) <- rownames(granu_de)

fgseaRes <- fgsea(pathways = msigdbr_list, 
                  stats = ranks)

result <- apply(fgseaRes,2,as.character)
write.csv(result, file = "OT3 Granulocyte CreERT2 VS CreERT2;Il33 fGSEA Result.csv")

fgseaRes_sig <- subset(fgseaRes, fgseaRes$padj <= 0.05)

waterfall_plot_sig <- function (fgseaRes, graph_title) {
  (fgseaRes %>% 
     mutate(short_name = str_split_fixed(pathway, "_",2)[,2])%>%
     ggplot( aes(reorder(short_name,NES), NES)) +
     geom_bar(stat= "identity", aes(fill = padj))+
     coord_flip()+
     labs(x = "Hallmark Gene Set", y = "Normalized Enrichment Score", title = graph_title)+
     theme(axis.text.y = element_text(size = 7), 
           plot.title = element_text(hjust = 1)))+scale_fill_viridis(direction = 1, option = "F")
}
waterfall_plot_sig(fgseaRes_sig, "Pathways Enriched in WT vs. KO - Granulocytes")


#### Figure S5D, T Cell Recruitment Heatmap ####
cd8_genes <- c("Cxcl9", "Ccl7", "Ccl2", 
               "Cxcr3", "Ccr1", "Ccr2", "Ccr3", "Ccr5", "Ccr10", "Ccr4")
avgexp <- AverageExpression(IL33_OT, return.seurat = T, features = cd8_genes, group.by = "Group")
levels(avgexp@active.ident)
DoHeatmap(avgexp, features = cd8_genes, slot = 'data', draw.lines = F) + scale_fill_viridis(direction = -1, option = "G")