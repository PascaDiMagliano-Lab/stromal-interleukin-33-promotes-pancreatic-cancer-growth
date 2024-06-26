#This script recreates the iKras scRNAseq analysis from Figures 5C-F and S6A-D of Donahue et al. 2024.
#This analysis utilizes datasets generated in this study, as well as previously published datasets from:

#iKrasG12D 3 weeks ON & iKrasG12D;Trp53R172H/+ orthotopic ON (NCBI GEO GSE140628)
#Zhang, Y., et al., Regulatory T-cell Depletion Alters the Tumor Microenvironment and Accelerates Pancreatic Carcinogenesis. 
#Cancer Discov, 2020. 10(3): p. 422-439.

#iKrasG12D 3 weeks OFF (NCBI GEO GSE179846)
#Velez-Delgado, A., et al., Extrinsic KRAS Signaling Shapes the Pancreatic Microenvironment Through Fibroblast Reprogramming. 
#Cell Mol Gastroenterol Hepatol, 2022. 13(6): p. 1673-1699.

#The remaining iKras datasets used in this analysis are reported for the first time in this study and are available at NCBI GEO GSE269888.

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

#### Preprocessing - Global Object ------------------------------------------------------------------------------------------------------------ ####
#Load Raw Data:
iKras_p53_on_14wk1_data <- Read10X_h5("~/iKras p53 14wk ON 1/raw_feature_bc_matrix.h5")
iKras_p53_on_14wk2_data <- Read10X_h5("~/iKras p53 14wk ON 2/raw_feature_bc_matrix.h5")
iKras_p53_on_15wk_off_1wk_data_1 <- Read10X_h5("~/15wk on 1wk OFF Nov 2019 1/raw_feature_bc_matrix.h5")
iKras_p53_on_15wk_off_1wk_data_2 <- Read10X_h5("~/15wk on 1wk OFF Nov 2019 2/raw_feature_bc_matrix.h5")

iKras_on_3wk_data <- Read10X("~/iKras 3 week ON/raw_feature_bc_matrix/")
iKras_on_3wk_off_3d_data <- Read10X("~/iKras 3 week OFF/raw_feature_bc_matrix/")

iKras_on_5wk_data <- Read10X_h5("~/iKras 5 week ON/raw_feature_bc_matrix.h5")
iKras_on_5wk_off_3d_data <- Read10X_h5("~/iKras 5 week OFF/raw_feature_bc_matrix.h5")

iKras_on_OT_data <- Read10X("~/OT 3 week ON/raw_feature_bc_matrix/")
iKras_on_OT_off_3d_data <- Read10X("~/OT 3 week OFF/raw_feature_bc_matrix/")

#Create Seurat Objects:
iKras_p53_14wk_1 <- CreateSeuratObject(iKras_p53_on_14wk1_data, min.cells = 3, min.features = 100)
iKras_p53_14wk_2 <- CreateSeuratObject(iKras_p53_on_14wk2_data, min.cells = 3, min.features = 100)
iKras_p53_15wk_1wk_off_1 <- CreateSeuratObject(iKras_p53_on_15wk_off_1wk_data_1, min.cells = 3, min.features = 100)
iKras_p53_15wk_1wk_off_2 <- CreateSeuratObject(iKras_p53_on_15wk_off_1wk_data_2, min.cells = 3, min.features = 100)

iKras_3wk <- CreateSeuratObject(iKras_on_3wk_data, min.cells = 3, min.features = 100)
iKras_3wk_off <- CreateSeuratObject(iKras_on_3wk_off_3d_data, min.cells = 3, min.features = 100)

iKras_5wk <- CreateSeuratObject(iKras_on_5wk_data, min.cells = 3, min.features = 100)
iKras_5wk_off <- CreateSeuratObject(iKras_on_5wk_off_3d_data, min.cells = 3, min.features = 100)

iKras_OT <- CreateSeuratObject(iKras_on_OT_data, min.cells = 3, min.features = 100)
iKras_OT_off <- CreateSeuratObject(iKras_on_OT_off_3d_data, min.cells = 3, min.features = 100)

#Add Desired Metadata:

iKras_p53_14wk_1[["KRAS_SIMPLE"]] <- "ON"
iKras_p53_14wk_2[["KRAS_SIMPLE"]] <- "ON"
iKras_p53_15wk_1wk_off_1[["KRAS_SIMPLE"]] <- "OFF"
iKras_p53_15wk_1wk_off_2[["KRAS_SIMPLE"]] <- "OFF"
iKras_3wk[["KRAS_SIMPLE"]] <- "ON"
iKras_3wk_off[["KRAS_SIMPLE"]] <- "OFF"
iKras_5wk[["KRAS_SIMPLE"]] <- "ON"
iKras_5wk_off[["KRAS_SIMPLE"]] <- "OFF"
iKras_OT[["KRAS_SIMPLE"]] <- "ON"
iKras_OT_off[["KRAS_SIMPLE"]] <- "OFF"

iKras_p53_14wk_1[["Time"]] <- "14wk_ON"
iKras_p53_14wk_2[["Time"]] <- "14wk_ON"
iKras_p53_15wk_1wk_off_1[["Time"]] <- "15wk_ON_1wk_OFF"
iKras_p53_15wk_1wk_off_2[["Time"]] <- "15wk_ON_1wk_OFF"
iKras_3wk[["Time"]] <- "3wk_ON"
iKras_3wk_off[["Time"]] <- "3wk_ON_3d_OFF"
iKras_5wk[["Time"]] <- "5wk_ON"
iKras_5wk_off[["Time"]] <- "5wk_ON_3d_OFF"
iKras_OT[["Time"]] <- "3wk_OT"
iKras_OT_off[["Time"]] <- "3wk_OT_3d_OFF"

iKras_p53_14wk_1[["Run_Date"]] <- "jun_2020"
iKras_p53_14wk_2[["Run_Date"]] <- "jun_2020"
iKras_p53_15wk_1wk_off_1[["Run_Date"]] <- "nov_2019"
iKras_p53_15wk_1wk_off_2[["Run_Date"]] <- "jun_2021"
iKras_3wk[["Run_Date"]] <- "mar_2019"
iKras_3wk_off[["Run_Date"]] <- "mar_2019"
iKras_5wk[["Run_Date"]] <- "oct_2019"
iKras_5wk_off[["Run_Date"]] <- "oct_2019"
iKras_OT[["Run_Date"]] <- "2018"
iKras_OT_off[["Run_Date"]] <- "2018"

#Merge Seurat Objects:
iKras <- merge(x = iKras_p53_14wk_1, y = c(iKras_p53_14wk_2,
                                                iKras_p53_15wk_1wk_off_1,
                                                iKras_p53_15wk_1wk_off_2,
                                                iKras_3wk,
                                                iKras_3wk_off,
                                                iKras_5wk,
                                                iKras_5wk_off,
                                                iKras_OT,
                                                iKras_OT_off), add.cell.ids = (c("a", 
                                                                                "b", 
                                                                                "c", 
                                                                                "d", 
                                                                                "e", 
                                                                                "f",
                                                                                "g",
                                                                                "h",
                                                                                "i",
                                                                                "j")))

#Normalize Data:
iKras <- NormalizeData(object = iKras, normalization.method = "LogNormalize", scale.factor = 10000)

#Apply Unbiased QC Cutoffs:
iKras[["percent.mt"]] <- PercentageFeatureSet(object = iKras, pattern = "^mt-")
iKras <- subset(x = iKras, subset = nCount_RNA > 800 & nCount_RNA < 100000 & percent.mt < 15)

#Integrate object by Run_Date:
Idents(object = iKras) <- "Run_Date"
iKras.list <- SplitObject(iKras, split.by = "Run_Date")

iKras.list <- lapply(X = iKras.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

iKras.anchors <- FindIntegrationAnchors(object.list = iKras.list)
iKras <- IntegrateData(anchorset = iKras.anchors)
DefaultAssay(iKras) <- "integrated"

#Scale Data:
iKras <- ScaleData(iKras, verbose = FALSE, features = rownames(iKras))

#Run PCA and Determine Dimensions for 90% Variance:
iKras <- RunPCA(iKras, npcs = 30, verbose = FALSE)
st_dev <- iKras@reductions$pca@stdev
var <- st_dev^2
sum(var[1:16])/ sum(var)

#Find Neighbors and Cluster Cells:
iKras <- FindNeighbors(object = iKras, dims = 1:16)
iKras <- FindClusters(object = iKras, resolution = 3)

#Run UMAP:
iKras <- RunUMAP(iKras, reduction = "pca", dims = 1:16, verbose = F)
DimPlot(iKras, reduction = "umap", label = TRUE)

#Identify Cell Populations:
DotPlot(iKras, features = c("Krt18", "Egfp","Msln","Lrrn4","Try4", "Col1a2", "Acta2", "Clec3b", "Cspg4","Ptprc", "Cd3e", "Trdc","Cd4", "Foxp3", "Ccr7", "Cd8a", "Nkg7","Il2ra","Il1rl1","Arg1", "Retnla","Cd68","Fcgr3","Adgre1", "Itgam","Cd14","Mrc1", "S100a8", "Cd33","H2-Eb1",
                                       "Itgae","Clec9a","Batf3", "Ly6c2", "Cd79a", "Cd19", "Ms4a1","Ccr10","Prdm1","Kit","Ccna2", "Pecam1", "Cdh5", "Hbb-bt"), cols = "RdYlBu", dot.scale = 8, assay = 'RNA') + RotatedAxis()

#Annotate Clusters:
Idents(iKras) <- "seurat_clusters"
iKras <- RenameIdents(iKras,
                                 "0" = "B Cell", 
                                 "1" = "iCAF",  
                                 "2" = "CD4+ T Cell", 
                                 "3" = "iCAF", 
                                 "4" = "iCAF", 
                                 "5" = "iCAF", 
                                 "6" = "CD8+ T Cell", 
                                 "7" = "CD8+ T Cell", 
                                 "8" = "CD4+ T Cell",
                                 "9" = "RBC", 
                                 "10" = "iCAF", 
                                 "11" = "iCAF", 
                                 "12" = "iCAF", 
                                 "13" = "CD4+ T Cell", 
                                 "14" = "iCAF", 
                                 "15" = "iCAF", 
                                 "16" = "iCAF", 
                                 "17" = "iCAF",
                                 "18" = "iCAF", 
                                 "19" = "RBC", 
                                 "20" = "iCAF", 
                                 "21" = "RBC",
                                 "22" = "Treg", 
                                 "23" = "Plasma Cell", 
                                 "24" = "iCAF", 
                                 "25" = "Macrophage", 
                                 "26" = "B Cell", 
                                 "27" = "RBC", 
                                 "28" = "Treg",
                                 "29" = "iCAF",
                                 "30" = "Granulocyte",
                                 "31" = "RBC",
                                 "32" = "Plasma Cell",
                                 "33" = "NK Cell",
                                 "34" = "yd T Cell", 
                                 "35" = "Proliferating Lymphocyte", 
                                 "36" = "myCAF",
                                 "37" = "iCAF",
                                 "38" = "NK T Cell",
                                 "39" = "Mesothelial",
                                 "40" = "iCAF",
                                 "41" = "Acinar",
                                 "42" = "CD8+ T Cell",
                                 "43" = "Mesothelial", 
                                 "44" = "Mast Cell",
                                 "45" = "ILC2", 
                                 "46" = "Granulocyte", 
                                 "47" = "cDC1",
                                 "48" = "RBC",
                                 "49" = "moDC",
                                 "50" = "Activated DC",
                                 "51" = "Endothelial",
                                 "52" = "iCAF",
                                 "53" = "Proliferating CAF",
                                 "54" = "pDC",
                                 "55" = "Plasma Cell",
                                 "56" = "iCAF",
                                 "57" = "Proliferating Myeloid",
                                 "58" = "RBC",
                                 "59" = "CK19+ ECAD+")
iKras[["detailed_clusters"]] <- iKras@active.ident

Idents(iKras) <- "detailed_clusters"
iKras <- RenameIdents(iKras,
                                 "iCAF" = "Fibroblast",
                                 "myCAF" = "Fibroblast",
                                 "Acinar" = "Acinar",
                                 "Mesothelial" = "Mesothelial",
                                 "Endothelial" = "Endothelial",
                                 "CK19+ ECAD+" = "CK19+ ECAD+",
                                 "CD4+ T Cell" = "CD4+ T Cell",
                                 "CD8+ T Cell" = "CD8+ T Cell",
                                 "Treg" = "CD4+ T Cell",
                                 "yd T Cell" = "yd T Cell",
                                 "ILC2" = "ILC2",
                                 "NK Cell" = "NK Cell",
                                 "NK T Cell" = "NK T Cell",
                                 "Macrophage" = "Macrophage",
                                 "Granulocyte" = "Granulocyte",
                                 "cDC1" = "Dendritic Cell",
                                 "moDC" = "Dendritic Cell",
                                 "Activated DC" = "Dendritic Cell",
                                 "pDC" = "Dendritic Cell",
                                 "Mast Cell" = "Mast Cell",
                                 "B Cell" = "B Cell", 
                                 "Plasma Cell" = "Plasma Cell", 
                                 "Proliferating Lymphocyte" = "Proliferating", 
                                 "Proliferating CAF" = "Proliferating",
                                 "Proliferating Myeloid" = "Proliferating",
                                 "RBC" = "RBC")

new_order <- c("CK19+ ECAD+",
               "Acinar",
               "Endothelial",
               "Mesothelial",
               "Fibroblast",
               "Macrophage",
               "Dendritic Cell",
               "Granulocyte",
               "NK Cell",
               "NK T Cell",
               "CD8+ T Cell",
               "ILC2",
               "yd T Cell",
               "CD4+ T Cell",
               "Mast Cell",
               "B Cell",
               "Plasma Cell",
               "Proliferating",
               "RBC")
iKras@active.ident <- factor(iKras@active.ident, levels = new_order)
iKras[["simple_clusters"]] <- iKras@active.ident

Idents(iKras) <- "simple_clusters"
iKras <- RenameIdents(iKras,
                      "CK19+ ECAD+" = "CK19+ ECAD+",
                      "Fibroblast" = "Fibroblast",
                      "Acinar" = "Acinar",
                      "Mesothelial" = "Mesothelial",
                      "Endothelial" = "Endothelial",
                      "Macrophage" = "Myeloid",
                      "Dendritic Cell" = "Myeloid",
                      "CD4+ T Cell" = "T Lymphocyte",
                      "CD8+ T Cell" = "T Lymphocyte",
                      "Granulocyte" = "Myeloid",
                      "yd T Cell" = "T Lymphocyte",
                      "ILC2" = "ILC2",
                      "NK Cell" = "NK Cell",
                      "NK T Cell" = "T Lymphocyte",
                      "Mast Cell" = "Myeloid",
                      "B Cell" = "B Lymphocyte", 
                      "Plasma Cell" = "B Lymphocyte", 
                      "Proliferating" = "Proliferating", 
                      "RBC" = "RBC")

new_order <- c("CK19+ ECAD+",
               "Acinar",
               "Endothelial",
               "Mesothelial",
               "Fibroblast",
               "Myeloid",
               "NK Cell",
               "T Lymphocyte",
               "ILC2",
               "B Lymphocyte",
               "Proliferating",
               "RBC")
iKras@active.ident <- factor(iKras@active.ident, levels = new_order)
iKras[["figure_clusters"]] <- iKras@active.ident

#Reorder Metadatas:
Idents(iKras) <- "KRAS_SIMPLE"
new_order <- c("ON",
               "OFF")
iKras@active.ident <- factor(iKras@active.ident, levels = new_order)
iKras[["KRAS_SIMPLE"]] <- iKras@active.ident

Idents(iKras) <- "Time"
new_order <- c("3wk_ON",
               "3wk_ON_3d_OFF",
               "5wk_ON",
               "5wk_ON_3d_OFF",
               "14wk_ON",
               "15wk_ON_1wk_OFF",
               "3wk_OT",
               "3wk_OT_3d_OFF")
iKras@active.ident <- factor(iKras@active.ident, levels = new_order)
iKras[["Time"]] <- iKras@active.ident

#Save Seurat Object:
Idents(iKras) <- "figure_clusters"
save(iKras, file = 'iKras_panin_and_OT.RData')

#### Preprocessing - Fibroblasts -------------------------------------------------------------------------------------------------------------- ####
#Subset and Re-Cluster Fibroblasts Alone:
Idents(iKras) <- "detailed_clusters"
ikras_fibros <- subset(iKras, idents = c("iCAF", "myCAF", "Proliferating CAF"))
DefaultAssay(ikras_fibros) <- "RNA"

Idents(object = ikras_fibros) <- "Run_Date"
table(ikras_fibros@active.ident)
#jun_2020 nov_2019 jun_2021 mar_2019 oct_2019     2018 
#    2474       83     4196     2378     3800     1561 

#Integrate object by Run_Date:
iKras.list <- SplitObject(ikras_fibros, split.by = "Run_Date")
iKras.list <- lapply(X = iKras.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

iKras.anchors <- FindIntegrationAnchors(object.list = iKras.list, k.filter = 83)
ikras_fibros <- IntegrateData(anchorset = iKras.anchors, k.weight = 83)
DefaultAssay(ikras_fibros) <- "integrated"

#Scale Data:
ikras_fibros <- ScaleData(ikras_fibros, verbose = FALSE, features = rownames(ikras_fibros))

#Run PCA and Determine Dimensions for 90% Variance:
ikras_fibros <- RunPCA(ikras_fibros, npcs = 30, verbose = FALSE)
st_dev <- ikras_fibros@reductions$pca@stdev
var <- st_dev^2
sum(var[1:21])/ sum(var)

#Find Neighbors and Cluster Cells:
ikras_fibros <- FindNeighbors(object = ikras_fibros, dims = 1:21)
ikras_fibros <- FindClusters(object = ikras_fibros, resolution = 3)

#Run UMAP:
ikras_fibros <- RunUMAP(ikras_fibros, reduction = "pca", dims = 1:21, verbose = F)
DimPlot(ikras_fibros, reduction = "umap", label = TRUE)

#Identify Cell Populations:
DotPlot(ikras_fibros, features = c("Krt18", "Egfp","Msln","Lrrn4","Try4", "Col1a2", "Acta2", "Clec3b", "Cspg4","Ptprc", "Cd3e", "Trdc","Cd4", "Foxp3", "Ccr7", "Cd8a", "Nkg7","Il2ra","Il1rl1","Arg1", "Retnla","Cd68","Fcgr3","Adgre1", "Itgam","Cd14","Mrc1", "S100a8", "Cd33","H2-Eb1",
                                   "Itgae","Clec9a","Batf3", "Ly6c2", "Cd79a", "Cd19", "Ms4a1","Ccr10","Prdm1","Kit","Ccna2", "Pecam1", "Cdh5", "Hbb-bt"), cols = "RdYlBu", dot.scale = 8, assay = 'RNA') + RotatedAxis()

#Remove Non-Fibroblast Cells and Repeat Integration:
ikras_fibros <- subset(ikras_fibros, idents = c("22", "31", "32", "33"), invert = T)

DefaultAssay(ikras_fibros) <- "RNA"
Idents(object = ikras_fibros) <- "Run_Date"
table(ikras_fibros@active.ident)
#jun_2020 nov_2019 jun_2021 mar_2019 oct_2019     2018 
#    2314       82     3950     2180     3647     1513   

iKras.list <- SplitObject(ikras_fibros, split.by = "Run_Date")
iKras.list <- lapply(X = iKras.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

iKras.anchors <- FindIntegrationAnchors(object.list = iKras.list, k.filter = 82)
ikras_fibros <- IntegrateData(anchorset = iKras.anchors, k.weight = 82)
DefaultAssay(ikras_fibros) <- "integrated"

#Scale Data:
ikras_fibros <- ScaleData(ikras_fibros, verbose = FALSE, features = rownames(ikras_fibros))

#Run PCA and Determine Dimensions for 90% Variance:
ikras_fibros <- RunPCA(ikras_fibros, npcs = 30, verbose = FALSE)
st_dev <- ikras_fibros@reductions$pca@stdev
var <- st_dev^2
sum(var[1:22])/ sum(var)

#Find Neighbors and Cluster Cells:
ikras_fibros <- FindNeighbors(object = ikras_fibros, dims = 1:22)
ikras_fibros <- FindClusters(object = ikras_fibros, resolution = 4)

#Run UMAP:
ikras_fibros <- RunUMAP(ikras_fibros, reduction = "pca", dims = 1:22, verbose = F)
DimPlot(ikras_fibros, reduction = "umap", label = TRUE)

#Ensure Purity of Fibroblasts:
DotPlot(ikras_fibros, features = c("Krt18", "Egfp","Msln","Lrrn4","Try4", "Col1a2", "Acta2", "Clec3b", "Cspg4","Ptprc", "Cd3e", "Trdc","Cd4", "Foxp3", "Ccr7", "Cd8a", "Nkg7","Il2ra","Il1rl1","Arg1", "Retnla","Cd68","Fcgr3","Adgre1", "Itgam","Cd14","Mrc1", "S100a8", "Cd33","H2-Eb1",
                                   "Itgae","Clec9a","Batf3", "Ly6c2", "Cd79a", "Cd19", "Ms4a1","Ccr10","Prdm1","Kit","Ccna2", "Pecam1", "Cdh5", "Hbb-bt"), cols = "RdYlBu", dot.scale = 8, assay = 'RNA') + RotatedAxis()

#Save Fibroblast Object:
save(ikras_fibros, file = 'fibros_ikras_panin_and_OT.RData')

#### Figures----------------------------------------------------------------------------------------------------------------------------------- ####
#Remove Proliferating CAFs for Differential Expression Analysis:
fibro_for_de <- subset(ikras_fibros, idents = "39", invert = T)

#### Figure 5C, Split Global UMAPs #### 
DimPlot(iKras, reduction = "umap", label = F, cols = c("CK19+ ECAD+" = "#b73a3a",
                                                       "Fibroblast" = "#80b1d3",
                                                       "Acinar" = "#108372",
                                                       "Mesothelial" = "#fee391",
                                                       "Endothelial" = "navy",
                                                       "T Lymphocyte" = "#8c6bb1",
                                                       "ILC2" = "#fccde5",
                                                       "NK Cell" = "#8dd3c7",
                                                       "Myeloid" = "#fdb462",
                                                       "B Lymphocyte" = "#bebada",
                                                       "Proliferating" = "#d9d9d9",
                                                       "RBC" = "#ccebc5"), pt.size = .5, group.by = "figure_clusters")

DimPlot(iKras, reduction = "umap", label = F, group.by = "Time",
        cols = c("5wk_ON" = "#648FFF",
                 "5wk_ON_3d_OFF" = "#648FFF", 
                 "3wk_ON" = "#DC267F",
                 "3wk_ON_3d_OFF" = "#DC267F",
                 "3wk_OT" = "#004D40",
                 "3wk_OT_3d_OFF" = "#004D40",
                 "14wk_ON" = "#FFB000",
                 "15wk_ON_1wk_OFF" = "#FFB000"), split.by = "KRAS_SIMPLE")

#### Figure 5D, Split Il33 Feature Plots #### 
DefaultAssay(iKras) <- "RNA"
FeaturePlot(iKras, features = "Il33", order = T, cols = c("gainsboro", "firebrick1"), split.by = "KRAS_SIMPLE", pt.size = 0.5)
DefaultAssay(iKras) <- "integrated"

#### Figure S6A, Global Cluster Marker Dotplot #### 
Idents(iKras) <- "figure_clusters"
DotPlot(iKras, features = c("Krt19", "Epcam", "Cdh1","Egfp",
                            "Try4", "Amy2a2",
                            "Pecam1", "Cdh5",
                            "Msln","Lrrn4",
                            "Col1a2", "Pdpn",
                            "Ptprc",
                            "Cd68","Adgre1", "Itgam","Cd14","H2-Eb1", "Batf3",
                            "Nkg7", "Gzmb",
                            "Cd3e", "Cd4", "Cd8a",
                            "Il2ra","Il1rl1","Arg1",
                            "Cd79a", "Cd19", "Ms4a1","Ccr10",
                            "Ccna2", "Mki67",
                            "Hbb-bt", "Hba-a1"), cols = "RdYlBu", dot.scale = 8, assay = 'RNA') + RotatedAxis()

#### Figure S6B, Cell Abundance Bar Graph #### 
Idents(object = iKras) <- 'Time'
table(iKras@active.ident)
# 3wk_ON   3wk_ON_3d_OFF          5wk_ON   5wk_ON_3d_OFF         14wk_ON 15wk_ON_1wk_OFF          3wk_OT 
#   2655            4333            4697            2017            5908            9015            1878 
#3wk_OT_3d_OFF 
#         2913

samples.list <- unique(iKras$Time)
clusters <- lapply(samples.list, function(x){
  subset <- subset(iKras, subset = Time == x)
  dist <- data.frame(table(subset$figure_clusters))
  
  return(dist)
})

names(clusters) <- samples.list

clusters_percent <- lapply(clusters, FUN = function(x){
  summ <- sum(x$Freq)
  x$Freq <- (x$Freq/summ)
  return(x)
})

clusters_dist <- reshape2::melt(clusters, id.var = "Var1")
colnames(clusters_dist) <- c("figure_clusters","variable","value","Time")
clusters_percent_dist <- reshape2::melt(clusters_percent, id.var = "Var1")
colnames(clusters_percent_dist) <- c("figure_clusters","variable","value","Time")

ggplot(clusters_percent_dist, aes(fill=figure_clusters, y = value, x = Time)) + 
  scale_x_discrete(limits = c("3wk_ON", "3wk_ON_3d_OFF", "5wk_ON", "5wk_ON_3d_OFF", "14wk_ON", "15wk_ON_1wk_OFF","3wk_OT", "3wk_OT_3d_OFF")) + 
  geom_bar(position = "stack", stat = "identity") + 
  theme_bw()+
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, size = 8, face = "bold", hjust = 1))+
  ylab("count") + xlab("orig.ident") + ggtitle("Relative Cell Types Abundance")+
  scale_fill_manual(values = c("CK19+ ECAD+" = "#b73a3a",
                               "Fibroblast" = "#80b1d3",
                               "Acinar" = "#108372",
                               "Mesothelial" = "#fee391",
                               "Endothelial" = "navy",
                               "T Lymphocyte" = "#8c6bb1",
                               "ILC2" = "#fccde5",
                               "NK Cell" = "#8dd3c7",
                               "Myeloid" = "#fdb462",
                               "B Lymphocyte" = "#bebada",
                               "Proliferating" = "#d9d9d9",
                               "RBC" = "#ccebc5"))

#### Figures 5E&S6C, Fibroblast Violin Plots ####
fibro_3_de <- FindMarkers(fibro_for_de, ident.1 = "3wk_ON", ident.2 = "3wk_ON_3d_OFF", group.by = "Time", assay = 'RNA')
fibro_3_de[["genes"]] <- row.names(fibro_3_de)
#p adj 1.311659e-18 Il33
#p adj 1.851406e-10 Igf1
#p adj 4.772939e-15 Timp2

fibro_5_de <- FindMarkers(fibro_for_de, ident.1 = "5wk_ON", ident.2 = "5wk_ON_3d_OFF", group.by = "Time", assay = 'RNA')
fibro_5_de[["genes"]] <- row.names(fibro_5_de)
#p adj 7.327054e-64 Il33
#p adj 3.200281e-170 Igf1
#p adj 3.873101e-07 Timp2

fibro_14_de <- FindMarkers(fibro_for_de, ident.1 = "14wk_ON", ident.2 = "15wk_ON_1wk_OFF", group.by = "Time", assay = 'RNA')
fibro_14_de[["genes"]] <- row.names(fibro_14_de)
#p adj 1.408223e-166 Il33
#p adj 1 Igf1
#p adj 6.793203e-40 Timp2

fibro_ot_de <- FindMarkers(fibro_for_de, ident.1 = "3wk_OT", ident.2 = "3wk_OT_3d_OFF", group.by = "Time", assay = 'RNA')
fibro_ot_de[["genes"]] <- row.names(fibro_ot_de)
#p adj 2.244752e-15 Il33
#p adj na Igf1
#p adj na Timp2

VlnPlot(fibro_for_de, features = "Il33", assay = 'RNA', group.by = "Time", pt.size = 0,
        cols = c("3wk_ON" = "#DC267F",
                 "3wk_ON_3d_OFF" = "#DCACC3",
                 "5wk_ON" = "#648FFF",
                 "5wk_ON_3d_OFF" = "#A4BBF7", 
                 "14wk_ON" = "#FFB000",
                 "15wk_ON_1wk_OFF" = "#E8D5AB",
                 "3wk_OT" = "#004D40",
                 "3wk_OT_3d_OFF" = "#8C9E9B")) +NoLegend()

VlnPlot(fibro_for_de, features = "Igf1", assay = 'RNA', group.by = "Time", pt.size = 0,
        cols = c("3wk_ON" = "#DC267F",
                 "3wk_ON_3d_OFF" = "#DCACC3",
                 "5wk_ON" = "#648FFF",
                 "5wk_ON_3d_OFF" = "#A4BBF7", 
                 "14wk_ON" = "#FFB000",
                 "15wk_ON_1wk_OFF" = "#E8D5AB",
                 "3wk_OT" = "#004D40",
                 "3wk_OT_3d_OFF" = "#8C9E9B")) +NoLegend()

VlnPlot(fibro_for_de, features = "Timp2", assay = 'RNA', group.by = "Time", pt.size = 0,
        cols = c("3wk_ON" = "#DC267F",
                 "3wk_ON_3d_OFF" = "#DCACC3",
                 "5wk_ON" = "#648FFF",
                 "5wk_ON_3d_OFF" = "#A4BBF7", 
                 "14wk_ON" = "#FFB000",
                 "15wk_ON_1wk_OFF" = "#E8D5AB",
                 "3wk_OT" = "#004D40",
                 "3wk_OT_3d_OFF" = "#8C9E9B")) +NoLegend()

#### Figures S6D, Mesothelial Violin Plots ####
Idents(iKras) <- "detailed_clusters"
mesothelial <- subset(iKras, idents = "Mesothelial")

VlnPlot(mesothelial, features = "Il33", assay = 'RNA', group.by = "Time", pt.size = 0,
        cols = c("3wk_ON" = "#DC267F",
                 "3wk_ON_3d_OFF" = "#DCACC3",
                 "5wk_ON" = "#648FFF",
                 "5wk_ON_3d_OFF" = "#A4BBF7", 
                 "14wk_ON" = "#FFB000",
                 "15wk_ON_1wk_OFF" = "#E8D5AB",
                 "3wk_OT" = "#004D40",
                 "3wk_OT_3d_OFF" = "#8C9E9B")) +NoLegend()

meso_3_de <- FindMarkers(mesothelial, ident.1 = "3wk_ON", ident.2 = "3wk_ON_3d_OFF", group.by = "Time", assay = 'RNA')
meso_3_de[["genes"]] <- row.names(meso_3_de)
#adj p 1

meso_5_de <- FindMarkers(mesothelial, ident.1 = "5wk_ON", ident.2 = "5wk_ON_3d_OFF", group.by = "Time", assay = 'RNA') 
meso_5_de[["genes"]] <- row.names(meso_5_de)
#adj p 1

meso_14_de <- FindMarkers(mesothelial, ident.1 = "14wk_ON", ident.2 = "15wk_ON_1wk_OFF", group.by = "Time", assay = 'RNA')
meso_14_de[["genes"]] <- row.names(meso_14_de)
#adj p 1

meso_ot_de <- FindMarkers(mesothelial, ident.1 = "3wk_OT", ident.2 = "3wk_OT_3d_OFF", group.by = "Time", assay = 'RNA')
meso_ot_de[["genes"]] <- row.names(meso_ot_de)
#adj p 1

#### Figures 5F, Fibroblast Pathway Analysis ####
h_gene_sets = msigdbr(species = "mouse", category = "H")
msigdbr_list = split(x = h_gene_sets$gene_symbol, f = h_gene_sets$gs_name)

fgsea_3_de <- FindMarkers(fibro_for_de, ident.1 = "3wk_ON", ident.2 = "3wk_ON_3d_OFF", group.by = "Time", assay = 'RNA', logfc.threshold = 0, min.pct = 0)
fgsea_5_de <- FindMarkers(fibro_for_de, ident.1 = "5wk_ON", ident.2 = "5wk_ON_3d_OFF", group.by = "Time", assay = 'RNA', logfc.threshold = 0, min.pct = 0) 
fgsea_14_de <- FindMarkers(fibro_for_de, ident.1 = "14wk_ON", ident.2 = "15wk_ON_1wk_OFF", group.by = "Time", assay = 'RNA', logfc.threshold = 0, min.pct = 0)
fgsea_ot_de <- FindMarkers(fibro_for_de, ident.1 = "3wk_OT", ident.2 = "3wk_OT_3d_OFF", group.by = "Time", assay = 'RNA', logfc.threshold = 0, min.pct = 0)

#Perform Enrichment Analysis:

### 3 week ###
ranks_3 <- fgsea_3_de$avg_log2FC
names(ranks_3) <- rownames(fgsea_3_de)
fgseaRes_3 <- fgsea(pathways = msigdbr_list, 
                    stats = ranks_3)
result_3 <- apply(fgseaRes_3,2,as.character)
write.csv(result_3, file = "iKras 3wk fibroblast ON VS OFF fGSEA Result.csv")
plotEnrichment(msigdbr_list[["HALLMARK_IL6_JAK_STAT3_SIGNALING"]], ranks_3, ticksSize = 0.5) + labs(title = "3 Week - ON VS OFF", subtitle="HALLMARK_IL6_JAK_STAT3_SIGNALING")

### 5 week ###
ranks_5 <- fgsea_5_de$avg_log2FC
names(ranks_5) <- rownames(fgsea_5_de)
fgseaRes_5 <- fgsea(pathways = msigdbr_list, 
                    stats = ranks_5)
result_5 <- apply(fgseaRes_5,2,as.character)
write.csv(result_5, file = "iKras 5wk iCAF ON VS OFF fGSEA Result.csv")
plotEnrichment(msigdbr_list[["HALLMARK_IL6_JAK_STAT3_SIGNALING"]], ranks_5, ticksSize = 0.5) + labs(title = "iCAF 5 Week - ON VS OFF", subtitle="HALLMARK_IL6_JAK_STAT3_SIGNALING")

### 14 week ###
ranks_14 <- fgsea_14_de$avg_log2FC
names(ranks_14) <- rownames(fgsea_14_de)
fgseaRes_14 <- fgsea(pathways = msigdbr_list, 
                     stats = ranks_14)
result_14 <- apply(fgseaRes_14,2,as.character)
write.csv(result_14, file = "iKras 14wk fibroblast ON VS OFF fGSEA Result.csv")
plotEnrichment(msigdbr_list[["HALLMARK_IL6_JAK_STAT3_SIGNALING"]], ranks_14, ticksSize = 0.5) + labs(title = "14+ Week - ON VS OFF", subtitle="HALLMARK_IL6_JAK_STAT3_SIGNALING")

### OT ###
ranks_ot <- fgsea_ot_de$avg_log2FC
names(ranks_ot) <- rownames(fgsea_ot_de)
fgseaRes_ot <- fgsea(pathways = msigdbr_list, 
                     stats = ranks_ot)
result_ot <- apply(fgseaRes_ot,2,as.character)
write.csv(result_ot, file = "iKras OT iCAF ON VS OFF fGSEA Result.csv")
plotEnrichment(msigdbr_list[["HALLMARK_IL6_JAK_STAT3_SIGNALING"]], ranks_ot, ticksSize = 0.5) + labs(title = "iCAF OT PDA - ON VS OFF", subtitle="HALLMARK_IL6_JAK_STAT3_SIGNALING")
