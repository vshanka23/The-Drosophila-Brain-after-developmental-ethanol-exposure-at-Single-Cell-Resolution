Sys.setenv(RETICULATE_PYTHON = "/usr/bin/python3")
library(cowplot)
library(ggplot2)
library(dplyr)
library(Seurat)
library(rlang)
library(TopKLists)
library(readxl)
library(gtools)
library(openxlsx)
library(readr)
options(future.globals.maxSize= 53687091200)

Female_Reference_R1_Data <- Read10X(data.dir="/home/vijay/Secretariat_data/Palmetto_sync/Projects/scRNASeq_ETOH/1_F_R_1/outs/filtered_feature_bc_matrix/")
Female_ETOH_R1_Data <- Read10X(data.dir="/home/vijay/Secretariat_data/Palmetto_sync/Projects/scRNASeq_ETOH/2_F_E_1/outs/filtered_feature_bc_matrix/")
Female_Reference_R2_Data <- Read10X(data.dir="/home/vijay/Secretariat_data/Palmetto_sync/Projects/scRNASeq_ETOH/3_F_R_2/outs/filtered_feature_bc_matrix/")
Female_ETOH_R2_Data <- Read10X(data.dir="/home/vijay/Secretariat_data/Palmetto_sync/Projects/scRNASeq_ETOH/4_F_E_2/outs/filtered_feature_bc_matrix/")
Male_Reference_R1_Data <- Read10X(data.dir="/home/vijay/Secretariat_data/Palmetto_sync/Projects/scRNASeq_ETOH/5_M_R_1/outs/filtered_feature_bc_matrix/")
Male_ETOH_R1_Data <- Read10X(data.dir="/home/vijay/Secretariat_data/Palmetto_sync/Projects/scRNASeq_ETOH/6_M_E_1/outs/filtered_feature_bc_matrix/")
Male_Reference_R2_Data <- Read10X(data.dir="/home/vijay/Secretariat_data/Palmetto_sync/Projects/scRNASeq_ETOH/7_M_R_2/outs/filtered_feature_bc_matrix/")
Male_ETOH_R2_Data <- Read10X(data.dir="/home/vijay/Secretariat_data/Palmetto_sync/Projects/scRNASeq_ETOH/8_M_E_2/outs/filtered_feature_bc_matrix/")

Female_Reference_R1 <- CreateSeuratObject(counts=Female_Reference_R1_Data,project="Female_Reference",min.cells = 5)
Female_Reference_R2 <- CreateSeuratObject(counts=Female_Reference_R2_Data,project="Female_Reference",min.cells = 5)
Female_ETOH_R1 <- CreateSeuratObject(counts=Female_ETOH_R1_Data,project="Female_ETOH",min.cells = 5)
Female_ETOH_R2 <- CreateSeuratObject(counts=Female_ETOH_R2_Data,project="Female_ETOH",min.cells = 5)
Male_Reference_R1 <- CreateSeuratObject(counts=Male_Reference_R1_Data,project="Male_Reference",min.cells = 5)
Male_Reference_R2 <- CreateSeuratObject(counts=Male_Reference_R2_Data,project="Male_Reference",min.cells = 5)
Male_ETOH_R1 <- CreateSeuratObject(counts=Male_ETOH_R1_Data,project="Male_ETOH",min.cells = 5)
Male_ETOH_R2 <- CreateSeuratObject(counts=Male_ETOH_R2_Data,project="Male_ETOH",min.cells = 5)

Female_Reference_R1$stim <- "Reference"
Female_Reference_R2$stim <- "Reference"
Female_ETOH_R1$stim <- "ETOH"
Female_ETOH_R2$stim <- "ETOH"
Male_Reference_R1$stim <- "Reference"
Male_Reference_R2$stim <- "Reference"
Male_ETOH_R1$stim <- "ETOH"
Male_ETOH_R2$stim <- "ETOH"

Female_Reference_R1$gender_stim <- "female_reference"
Female_Reference_R2$gender_stim <- "female_reference"
Female_ETOH_R1$gender_stim <- "female_ETOH"
Female_ETOH_R2$gender_stim <- "female_ETOH"
Male_Reference_R1$gender_stim <- "male_reference"
Male_Reference_R2$gender_stim <- "male_reference"
Male_ETOH_R1$gender_stim <- "male_ETOH"
Male_ETOH_R2$gender_stim <- "male_ETOH"

Female_Reference_R1$sample_id <- "female_reference_R1"
Female_Reference_R2$sample_id <- "female_reference_R2"
Female_ETOH_R1$sample_id <- "female_ETOH_R1"
Female_ETOH_R2$sample_id <- "female_ETOH_R2"
Male_Reference_R1$sample_id <- "male_reference_R1"
Male_Reference_R2$sample_id <- "male_reference_R2"
Male_ETOH_R1$sample_id <- "male_ETOH_R1"
Male_ETOH_R2$sample_id <- "male_ETOH_R2"

VlnPlot(Female_Reference_R1,features=c("nFeature_RNA","nCount_RNA"), ncol=2)
VlnPlot(Female_Reference_R2,features=c("nFeature_RNA","nCount_RNA"), ncol=2)
VlnPlot(Female_ETOH_R1,features=c("nFeature_RNA","nCount_RNA"), ncol=2)
VlnPlot(Female_ETOH_R2,features=c("nFeature_RNA","nCount_RNA"), ncol=2)
VlnPlot(Male_Reference_R1,features=c("nFeature_RNA","nCount_RNA"), ncol=2)
VlnPlot(Male_Reference_R2,features=c("nFeature_RNA","nCount_RNA"), ncol=2)
VlnPlot(Male_ETOH_R1,features=c("nFeature_RNA","nCount_RNA"), ncol=2)
VlnPlot(Male_ETOH_R2,features=c("nFeature_RNA","nCount_RNA"), ncol=2)

Female_Reference_R1 <- subset(Female_Reference_R1,subset=nFeature_RNA > 300 & nFeature_RNA < 2500)
VlnPlot(Female_Reference_R1,features=c("nFeature_RNA","nCount_RNA"), ncol=2)
Female_Reference_R2 <- subset(Female_Reference_R2,subset=nFeature_RNA > 300 & nFeature_RNA < 2500)
VlnPlot(Female_Reference_R2,features=c("nFeature_RNA","nCount_RNA"), ncol=2)
Female_ETOH_R1 <- subset(Female_ETOH_R1,subset=nFeature_RNA > 300 & nFeature_RNA < 2500)
VlnPlot(Female_ETOH_R1,features=c("nFeature_RNA","nCount_RNA"), ncol=2)
Female_ETOH_R2 <- subset(Female_ETOH_R2,subset=nFeature_RNA > 300 & nFeature_RNA < 2500)
VlnPlot(Female_ETOH_R2,features=c("nFeature_RNA","nCount_RNA"), ncol=2)
Male_Reference_R1 <- subset(Male_Reference_R1,subset=nFeature_RNA > 300 & nFeature_RNA < 2500)
VlnPlot(Male_Reference_R1,features=c("nFeature_RNA","nCount_RNA"), ncol=2)
Male_Reference_R2 <- subset(Male_Reference_R2,subset=nFeature_RNA > 300 & nFeature_RNA < 2500)
VlnPlot(Male_Reference_R2,features=c("nFeature_RNA","nCount_RNA"), ncol=2)
Male_ETOH_R1 <- subset(Male_ETOH_R1,subset=nFeature_RNA > 300 & nFeature_RNA < 2500)
VlnPlot(Male_ETOH_R1,features=c("nFeature_RNA","nCount_RNA"), ncol=2)
Male_ETOH_R2 <- subset(Male_ETOH_R2,subset=nFeature_RNA > 300 & nFeature_RNA < 2500)
VlnPlot(Male_ETOH_R2,features=c("nFeature_RNA","nCount_RNA"), ncol=2)

Female_Reference_R1 <- SCTransform(Female_Reference_R1,verbose = FALSE,return.only.var.genes = FALSE)
Female_Reference_R2 <- SCTransform(Female_Reference_R2,verbose = FALSE,return.only.var.genes = FALSE)
Female_ETOH_R1 <- SCTransform(Female_ETOH_R1,verbose = FALSE,return.only.var.genes = FALSE)
Female_ETOH_R2 <- SCTransform(Female_ETOH_R2,verbose = FALSE,return.only.var.genes = FALSE)
Male_Reference_R1 <- SCTransform(Male_Reference_R1,verbose = FALSE,return.only.var.genes = FALSE)
Male_Reference_R2 <- SCTransform(Male_Reference_R2,verbose = FALSE,return.only.var.genes = FALSE)
Male_ETOH_R1 <- SCTransform(Male_ETOH_R1,verbose = FALSE,return.only.var.genes = FALSE)
Male_ETOH_R2 <- SCTransform(Male_ETOH_R2,verbose = FALSE,return.only.var.genes = FALSE)

Integration_feature_set <- SelectIntegrationFeatures(object.list = c(Female_Reference_R1,Female_Reference_R2,Female_ETOH_R1,Female_ETOH_R2,Male_Reference_R1,Male_Reference_R2,Male_ETOH_R1,Male_ETOH_R2),nfeatures = 1500)

Integration_list <- PrepSCTIntegration(object.list = c(Female_Reference_R1,Female_Reference_R2,Female_ETOH_R1,Female_ETOH_R2,Male_Reference_R1,Male_Reference_R2,Male_ETOH_R1,Male_ETOH_R2), anchor.features = Integration_feature_set, verbose=FALSE)

Integration_anchors <- FindIntegrationAnchors(object.list = Integration_list, normalization.method = "SCT", anchor.features = Integration_feature_set, verbose=FALSE)

Integrated <- IntegrateData(anchorset = Integration_anchors, normalization.method = "SCT", verbose=FALSE)
Integrated <- RunPCA(object = Integrated, verbose=FALSE)

ElbowPlot(Integrated)

Integrated <- RunUMAP(Integrated, reduction = "pca", dims= 1:10)

Integrated <- FindNeighbors(Integrated, reduction = "pca", dims = 1:10)

cluster <- FindClusters(Integrated, resolution = 0.4)
cluster <- FindClusters(Integrated, resolution = 0.5)
cluster <- FindClusters(Integrated, resolution = 0.6)
cluster <- FindClusters(Integrated, resolution = 0.7)
cluster <- FindClusters(Integrated, resolution = 0.8)
cluster <- FindClusters(Integrated, resolution = 0.9)
cluster <- FindClusters(Integrated, resolution = 1.0)
cluster <- FindClusters(Integrated, resolution = 1.1)
cluster <- FindClusters(Integrated, resolution = 1.2)
cluster <- FindClusters(Integrated, resolution = 1.3)
cluster <- FindClusters(Integrated, resolution = 1.4)
cluster <- FindClusters(Integrated, resolution = 1.5)
cluster <- FindClusters(Integrated, resolution = 1.6)
cluster <- FindClusters(Integrated, resolution = 1.7)
cluster <- FindClusters(Integrated, resolution = 1.8)
cluster <- FindClusters(Integrated, resolution = 1.9)
cluster <- FindClusters(Integrated, resolution = 2.0)

Integrated <- FindClusters(Integrated, resolution =0.9)

Integrated$celltype <- Idents(Integrated)
plots <- DimPlot(Integrated, reduction= "umap", group.by = c("stim","gender_stim","celltype","sample_id"),combine=FALSE)

DimPlot(Integrated, group.by = "sample_id", pt.size = 0.001)
DimPlot(Integrated, group.by = "gender_stim", pt.size = 0.001)
DimPlot(Integrated, group.by = "stim", pt.size = 0.001)
DimPlot(Integrated, group.by = "celltype", pt.size = 0.001)

Integrated$celltype.stim <- paste(Idents(Integrated),Integrated$stim,sep = "_")
Integrated$celltype.gender_stim <- paste(Idents(Integrated),Integrated$gender_stim,sep = "_")

DefaultAssay(Integrated) <- "RNA"
Integrated <- NormalizeData(Integrated,verbose=FALSE)
Integrated_markers_all <- FindAllMarkers(Integrated,min.pct = 0.25, logfc.threshold = 0.5, only.pos = TRUE)
Integrated_cluster_markers <- Integrated_markers_all %>% group_by(cluster) %>% top_n(n=3, wt = avg_logFC) %>% print(n=3*43)
View(Integrated_cluster_markers)
write.csv(Integrated_cluster_markers,"scETOH_cluster_markers_top3.csv")

DefaultAssay(Integrated) <- "integrated"
Idents(Integrated)<-"celltype.gender_stim"

gender <- c("female","male")
group <- c("reference","ETOH")
cluster <- c(0:42)

female_ETOH_index <- as.list(paste(cluster,rep(gender[1],43),rep(group[2],43),sep="_"))
female_reference_index <- as.list(paste(cluster,rep(gender[1],43),rep(group[1],43),sep="_"))
male_ETOH_index <- as.list(paste(cluster,rep(gender[2],43),rep(group[2],43),sep="_"))
male_reference_index <- as.list(paste(cluster,rep(gender[2],43),rep(group[1],43),sep="_"))
##Check one cluster's worth of DE just to confirm that the for loop works correctly.

for (i in 0:42) {
Female_DE <- paste("Female_DE","C",i,sep = "")
assign(Female_DE,FindMarkers(Integrated,ident.1 = female_ETOH_index[i+1],ident.2 = female_reference_index[i+1],test.use = "MAST",assay = "SCT",slot = "scale.data"))
}

for (i in 0:42) {
Male_DE <- paste("Male_DE","C",i,sep = "")
assign(Male_DE,FindMarkers(Integrated,ident.1 = male_ETOH_index[i+1],ident.2 = male_reference_index[i+1],test.use = "MAST",assay = "SCT",slot = "scale.data"))
}

fbgn_annotation_ID_fb_2020_06 <- read_delim("~/Downloads/fbgn_annotation_ID_fb_2020_06.tsv","\t", escape_double = FALSE, trim_ws = TRUE, skip = 3)

Female_DE_list <- mget(mixedsort(ls(pattern="Female_DEC")))
Female_DE_list_FBGN_symbol <- lapply(Female_DE_list, function(df) {
df$FBGN <- with(fbgn_annotation_ID_fb_2020_06, `primary_FBgn#`[match(sapply(strsplit(rownames(df),"-"),"[",2),`annotation_ID`)])
df$symbol <- with(fbgn_annotation_ID_fb_2020_06, `##gene_symbol`[match(sapply(strsplit(rownames(df),"-"),"[",2),`annotation_ID`)])
df
})

Male_DE_list <- mget(mixedsort(ls(pattern="Male_DEC")))
Male_DE_list_FBGN_symbol <- lapply(Male_DE_list, function(df) {
df$FBGN <- with(fbgn_annotation_ID_fb_2020_06, `primary_FBgn#`[match(sapply(strsplit(rownames(df),"-"),"[",2),`annotation_ID`)])
df$symbol <- with(fbgn_annotation_ID_fb_2020_06, `##gene_symbol`[match(sapply(strsplit(rownames(df),"-"),"[",2),`annotation_ID`)])
df
})

write.xlsx(Female_DE_list_FBGN_symbol,file="scETOH_female_DE.xlsx",rowNames=TRUE)
write.xlsx(Male_DE_list_FBGN_symbol,file="scETOH_male_DE.xlsx",rowNames=TRUE)
