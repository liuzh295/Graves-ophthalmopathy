library(dplyr)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(viridis)
library(harmony)
library(ArchR)
setwd("/data/TAO")
export PATH=/home/cellranger-7.1.0:$PATH
indir=/data/data1/TAO/data
outdir=/data/data3/TAO
threads=10
refpath=/data/data1/cellranger_scRNA_reference/refdata-gex-GRCh38-2020-A 
cellranger=/home/cellranger-7.1.0/cellranger
cd ${outdir}

sample=SCS-20  
${cellranger} count --id=${sample} \
--transcriptome=${refpath} \
--fastqs=${indir} \
--sample=${sample} \
--localcores=${threads}
sample=SCS-15  
${cellranger} count --id=${sample} \
--transcriptome=${refpath} \
--fastqs=${indir} \
--sample=${sample} \
--localcores=${threads}
sample=SCS-11 
${cellranger} count --id=${sample} \
--transcriptome=${refpath} \
--fastqs=${indir} \
--sample=${sample} \
--localcores=${threads}
sample=SCS-09  
${cellranger} count --id=${sample} \
--transcriptome=${refpath} \
--fastqs=${indir} \
--sample=${sample} \
--localcores=${threads}
sample=SCS-07 
${cellranger} count --id=${sample} \
--transcriptome=${refpath} \
--fastqs=${indir} \
--sample=${sample} \
--localcores=${threads}
sample=SCS-03 
${cellranger} count --id=${sample} \
--transcriptome=${refpath} \
--fastqs=${indir} \
--sample=${sample} \
--localcores=${threads}
sample=SCS-01  
${cellranger} count --id=${sample} \
--transcriptome=${refpath} \
--fastqs=${indir} \
--sample=${sample} \
--localcores=${threads}

sample=W-1 
${cellranger} count --id=${sample} \
--transcriptome=${refpath} \
--fastqs=${indir} \
--sample=${sample} \
--localcores=${threads}
sample=W-02 
${cellranger} count --id=${sample} \
--transcriptome=${refpath} \
--fastqs=${indir} \
--sample=${sample} \
--localcores=${threads}
sample=W-03 
${cellranger} count --id=${sample} \
--transcriptome=${refpath} \
--fastqs=${indir} \
--sample=${sample} \
--localcores=${threads}
sample=W-4 
${cellranger} count --id=${sample} \
--transcriptome=${refpath} \
--fastqs=${indir} \
--sample=${sample} \
--localcores=${threads}
sample=W-5 
${cellranger} count --id=${sample} \
--transcriptome=${refpath} \
--fastqs=${indir} \
--sample=${sample} \
--localcores=${threads}

#data processing
E_H1 <- Read10X(data.dir = "/data/TAO/SCS-20/outs/filtered_feature_bc_matrix")
E_H1 <- CreateSeuratObject(counts = E_H1, project = "E_H-01", min.cells = 3, min.features = 400)
E_H1[["percent.mt"]] <- PercentageFeatureSet(E_H1, pattern = "^MT-") 

E_H2 <- Read10X(data.dir = "/data/TAO/SCS-15/outs/filtered_feature_bc_matrix")
E_H2 <- CreateSeuratObject(counts = E_H2, project = "E_H-02", min.cells = 3, min.features = 400)
E_H2[["percent.mt"]] <- PercentageFeatureSet(E_H2, pattern = "^MT-") 

E_T1 <- Read10X(data.dir = "/data/TAO/SCS-1/outs/filtered_feature_bc_matrix")
E_T1 <- CreateSeuratObject(counts = E_T1, project = "E_T-01", min.cells = 3, min.features = 400)
E_T1[["percent.mt"]] <- PercentageFeatureSet(E_T1, pattern = "^MT-") 

E_T2 <- Read10X(data.dir = "/data/TAO/SCS-11/outs/filtered_feature_bc_matrix")
E_T2 <- CreateSeuratObject(counts = E_T2, project = "E_T-02", min.cells = 3, min.features = 400)
E_T2[["percent.mt"]] <- PercentageFeatureSet(E_T2, pattern = "^MT-") 

E_T3 <- Read10X(data.dir = "/data/TAO/SCS-9/outs/filtered_feature_bc_matrix")
E_T3 <- CreateSeuratObject(counts = E_T3, project = "E_T-03", min.cells = 3, min.features = 400)
E_T3[["percent.mt"]] <- PercentageFeatureSet(E_T3, pattern = "^MT-") 

E_T4 <- Read10X(data.dir = "/data/TAO/SCS-7/outs/filtered_feature_bc_matrix")
E_T4 <- CreateSeuratObject(counts = E_T4, project = "E_T-04", min.cells = 3, min.features = 400)
E_T4[["percent.mt"]] <- PercentageFeatureSet(E_T4, pattern = "^MT-") 

E_T5 <- Read10X(data.dir = "/data/TAO/SCS-3/outs/filtered_feature_bc_matrix")
E_T5 <- CreateSeuratObject(counts = E_T5, project = "E_T-05", min.cells = 3, min.features = 400)
E_T5[["percent.mt"]] <- PercentageFeatureSet(E_T5, pattern = "^MT-") 

E_T6 <- Read10X(data.dir = "/data/TAO/SCS-3/outs/filtered_feature_bc_matrix")
E_T6 <- CreateSeuratObject(counts = E_T6, project = "E_T-05", min.cells = 3, min.features = 400)
E_T6[["percent.mt"]] <- PercentageFeatureSet(E_T6, pattern = "^MT-") 

E_T7 <- Read10X(data.dir = "/data/TAO/SCS-3/outs/filtered_feature_bc_matrix")
E_T7 <- CreateSeuratObject(counts = E_T7, project = "E_T-05", min.cells = 3, min.features = 400)
E_T7[["percent.mt"]] <- PercentageFeatureSet(E_T7, pattern = "^MT-") 

E_T8 <- Read10X(data.dir = "/data/TAO/SCS-3/outs/filtered_feature_bc_matrix")
E_T8 <- CreateSeuratObject(counts = E_T8, project = "E_T-05", min.cells = 3, min.features = 400)
E_T8[["percent.mt"]] <- PercentageFeatureSet(E_T8, pattern = "^MT-") 

E_T9 <- Read10X(data.dir = "/data/TAO/SCS-3/outs/filtered_feature_bc_matrix")
E_T9 <- CreateSeuratObject(counts = E_T9, project = "E_T-05", min.cells = 3, min.features = 400)
E_T9[["percent.mt"]] <- PercentageFeatureSet(E_T9, pattern = "^MT-") 

E_T10 <- Read10X(data.dir = "/data/TAO/SCS-3/outs/filtered_feature_bc_matrix")
E_T10 <- CreateSeuratObject(counts = E_T10, project = "E_T-05", min.cells = 3, min.features = 400)
E_T10[["percent.mt"]] <- PercentageFeatureSet(E_T10, pattern = "^MT-") 

E_H1 <- subset(E_H1, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 20)
E_H2 <- subset(E_H2, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 20)

E_T1 <- subset(E_T1, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 20)
E_T2 <- subset(E_T2, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 20)
E_T3 <- subset(E_T3, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 20)
E_T4 <- subset(E_T4, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 20)
E_T5 <- subset(E_T5, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 20)
E_T6 <- subset(E_T6, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 20)
E_T7 <- subset(E_T7, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 20)
E_T8 <- subset(E_T8, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 20)
E_T9 <- subset(E_T9, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 20)
E_T10 <- subset(E_T10, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 20)


scRNA.integrated <- merge(E_H1, y=c(E_H2,E_T1,E_T2,E_T3,E_T4,E_T5,E_T6,E_T7,E_T8,E_T9,E_T10))
DefaultAssay(scRNA.integrated) <- "RNA"
scRNA.integrated <- NormalizeData(scRNA.integrated, normalization.method = "LogNormalize", scale.factor = 10000)
scRNA.integrated <- FindVariableFeatures(scRNA.integrated, nfeatures = 4000)
scRNA.integrated <- ScaleData(scRNA.integrated, verbose = FALSE)
scRNA.integrated <- RunPCA(scRNA.integrated, verbose = FALSE)
scRNA.integrated =  scRNA.integrated %>% RunHarmony("orig.ident", plot_convergence = TRUE)
scRNA.integrated <- RunUMAP(scRNA.integrated, reduction = "harmony", dims = 1:30)
scRNA.integrated <- FindNeighbors(scRNA.integrated, reduction = "harmony", dims = 1:30)

scRNA.integrated <- FindClusters(scRNA.integrated, resolution = 0.6)

Idents(scRNA.integrated) <- "orig.ident"
DiseaseStates <- c("E_H-01"="F_Healthy", 
                   "E_H-02"="F_Healthy", 
                   "E_T-01"="GO", 
                   "E_T-02"="GO",
                   "E_T-03"="GO", 
                   "E_T-04"="GO", 
                   "E_T-05"="GO"
)
scRNA.integrated <- RenameIdents(scRNA.integrated, DiseaseStates) 
scRNA.integrated$DiseaseStates <- scRNA.integrated@active.ident
Idents(scRNA.integrated) <- "DiseaseStates"
table(scRNA.integrated@meta.data$DiseaseStates)

new.cluster.ids <- c("LPF1",   
                     "MYF1",
                     "COF",
                     "LPF2",
                     "Immune cell",
                     "Pericytes",
                     "Endothelial",
                     "Immune cell",
                     "Immune cell",
                     "Immune cell",  
                     "LPF3", 
                     "Immune cell",
                     "Progenitor",
                     "Immune cell",
                     "Unknown",
                     "MYF2",
                     "Immune cell")
identities <- as.character(scRNA.integrated@meta.data$seurat_clusters)
for (i in 0:length(new.cluster.ids)){
  identities[identities==as.character(i)] <- new.cluster.ids[i+1]
}
scRNA.integrated <- AddMetaData(scRNA.integrated, identities, col.name = "celltype")
scRNA.integrated2 <- subset(scRNA.integrated, idents = c("Unknown","Immune cell"),invert=TRUE)
pdf(paste0("Disease_state1.pdf"), width = 6, onefile=F)
print(DimPlot(scRNA.integrated, reduction = "umap", group.by = "DiseaseStates"))
dev.off()

#
scRNA.integrated1 <- subset(scRNA.integrated, idents = c("Immune cell"))
DefaultAssay(scRNA.integrated1) <- "RNA"
scRNA.integrated1 <- NormalizeData(scRNA.integrated1, normalization.method = "LogNormalize", scale.factor = 10000)
scRNA.integrated1 <- FindVariableFeatures(scRNA.integrated1, nfeatures = 4000)
scRNA.integrated1 <- ScaleData(scRNA.integrated1,vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = FALSE)
scRNA.integrated1 <- RunPCA(scRNA.integrated1, verbose = FALSE)
scRNA.integrated1 <- scRNA.integrated1 %>% RunHarmony("orig.ident", plot_convergence = TRUE)
scRNA.integrated1 <- RunUMAP(scRNA.integrated1,reduction = "harmony", dims = 1:30)
scRNA.integrated1 <- FindNeighbors(scRNA.integrated1,reduction = "harmony", dims = 1:30)
new.cluster.ids <- c("TC",   
                     "TC",
                     "Monocyte",
                     "NK",
                     "TC",
                     "NK",
                     "BC",
                     "DC",
                     "BC",
                     "DC",  
                     "Monocyte",
                     "Progenitor",
                     "Cycling",
                     "Progenitor",
                     "Monocyte")
identities <- as.character(scRNA.integrated1@meta.data$seurat_clusters)
for (i in 0:length(new.cluster.ids)){
  identities[identities==as.character(i)] <- new.cluster.ids[i+1]
}
scRNA.integrated1 <- AddMetaData(scRNA.integrated1, identities, col.name = "celltype")
Idents(scRNA.integrated1) <- "celltype"
DefaultAssay(scRNA.integrated1) <- "RNA"
genes_to_check = c(
  'CD3D',   # T cell
  "CD4","CD8A",
  "GPR183", 
  "FOXP3", 
  "CCR7","SELL",   
  "NKG7","GNLY",  
  "GZMB","GZMA","PRF1",
  "MKI67",
  "TRGV9","TRDV2",       
  "TRAV1-2"          
)
p <- DotPlot(scRNA.integrated1, 
             features = genes_to_check, 
             assay = 'RNA',
             cols = c('lightgrey', '#8B0000')) +RotatedAxis()
ggsave("Immune_All_marker.pdf",p, width = 13, height = 8)
pdf(paste0("Disease_state.pdf"), width = 6, onefile=F)
print(DimPlot(scRNA.integrated1, reduction = "umap", group.by = "DiseaseState"))
dev.off()

#
scRNA.integrated <- merge(scRNA.integrated1, y=c(TC))
DefaultAssay(scRNA.integrated) <- "RNA"
scRNA.integrated <- NormalizeData(scRNA.integrated, normalization.method = "LogNormalize", scale.factor = 10000)
scRNA.integrated <- FindVariableFeatures(scRNA.integrated, nfeatures = 4000)
scRNA.integrated <- ScaleData(scRNA.integrated, vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = FALSE)
scRNA.integrated <- RunPCA(scRNA.integrated, npcs = 30, verbose = FALSE)
scRNA.integrated <- scRNA.integrated %>% RunHarmony("orig.ident", plot_convergence = TRUE)
scRNA.integrated <- RunUMAP(scRNA.integrated,reduction = "harmony", dims = 1:30)
scRNA.integrated <- FindNeighbors(scRNA.integrated,reduction = "harmony", dims = 1:30)
scRNA.integrated <- FindClusters(scRNA.integrated, resolution = 0.6)

Idents(scRNA.integrated) <- "seurat_clusters"
celltype <- c("0"="TC", 
              "1"="TC", 
              "2"="Monocyte", 
              "3"="TC",
              "4"="NK", 
              "5"="NK", 
              "6"="TC",
              "7"="BC",
              "8"="Unkown",
              "9"="TC", 
              "10"="BC", 
              "11"="TC", 
              "12"="Monocyte",
              "13"="TC", 
              "14"="TC", 
              "15"="TC",
              "16"="Monocyte",
              "17"="TC", 
              "11"="TC", 
              "18"="Monocyte",
              "19"="TC", 
              "20"="DC", 
              "21"="Unkown",
              "22"="Unkown",
              "23"="NK", 
              "24"="BC", 
              "25"="Monocyte",
              "26"="TC", 
              "27"="DC", 
              "28"="Cycling",
              "29"="BC",
              "30"="Monocyte",
              "31"="Unkown",
              "32"="Unkown",
              "33"="Progenitor"
)
scRNA.integrated <- RenameIdents(scRNA.integrated, celltype) 
scRNA.integrated$celltype <- scRNA.integrated@active.ident
Idents(scRNA.integrated) <- "celltype"
scRNA.integrated1 <- subset(scRNA.integrated, idents = c("Unkown"),invert=TRUE)

P<-VlnPlot(scRNA.integrated, features = c("IDNK","SLC35G1"),cols=cols,pt.size = 0.02)
ggsave("Immune_IDNK,SLC35G1.pdf",p, width = 13, height = 8)

pdf(paste0("Immune_Disease_state.pdf"), width = 6, onefile=F)
print(DimPlot(scRNA.integrated, reduction = "umap", group.by = "DiseaseState"))
dev.off()

Idents(scRNA.integrated) <- "celltype"
DefaultAssay(scRNA.integrated1) <- "RNA"
genes_to_check = c("CD3D","CD4", "CD8A", 
                   "NKG7","NCR1",  #NK
                   "CD79A","MS4A1",    #B cell
                   "CD14","S100A9", #Mo
                   "CD1C","CLEC4C", #DC 
                   "GATA2", #progenitor
                   "MKI67" #Cycling
)
P <- DotPlot(scRNA.integrated, 
             features = genes_to_check, 
             assay = 'RNA',
             cols = c('lightgrey', '#8B0000')) +RotatedAxis()
ggsave("Immune_RNA_marker.pdf",P, width = 7, height = 3)


library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
data.input = Tissue$data # normalized data matrix
meta = Tissue$meta # a dataframe with rownames containing cell mata data
cell.use = rownames(meta)[meta$condition == "LS"] # extract the cell names from disease data
data.input = data.input[, cell.use]
meta = meta[cell.use, ]
# meta = data.frame(labels = meta$labels[cell.use], row.names = colnames(data.input)) # manually create a dataframe consisting of the cell labels
unique(meta$labels) # check the cell labels
data.input = data.input[, cell.use]
meta = meta[cell.use, ]
# meta = data.frame(labels = meta$labels[cell.use], row.names = colnames(data.input)) # manually create a dataframe consisting of the cell labels
unique(meta$labels) # check the cell labels
data.input = data.input[, cell.use]
meta = meta[cell.use, ]
# meta = data.frame(labels = meta$labels[cell.use], row.names = colnames(data.input)) # manually create a dataframe consisting of the cell labels
unique(meta$labels) # check the cell labels
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
future::plan("multiprocess", workers = 4) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
levels(cellchat@idents)
p = netVisual_bubble(cellchat, sources.use = c(3,5,7,8,9), 
                     targets.use = c(1,2,4,6), remove.isolate = FALSE)
ggsave("bubble.pdf", p, width = 8, height = 12)

mat <- cellchat@net$count
par(mfrow = c(3,3), xpd=TRUE)
for (i in 1:nrow(mat)) {
  # i = 1
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, arrow.width = 0.2,
                   arrow.size = 0.1, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")

