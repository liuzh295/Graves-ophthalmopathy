library(dplyr)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(viridis)
library(ArchR)
library(harmony)
library(ggpubr)

setwd("~/scRNA")
scRNA <- readRDS("scRNA.rds")
DefaultAssay(scRNA) <- "RNA"
scRNA <- NormalizeData(scRNA, normalization.method = "LogNormalize", scale.factor= 1e4)
scRNA <- FindVariableFeatures(scRNA, nfeatures = 4000)                      
scRNA <- ScaleData(scRNA, vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = FALSE)
scRNA <- RunPCA(scRNA, npcs = 50, verbose = FALSE)
scRNA <- scRNA %>% RunHarmony("orig.ident", plot_convergence = TRUE)
scRNA <- RunUMAP(scRNA, reduction = "harmony", dims = 1:30)
scRNA <- FindNeighbors(scRNA, reduction = "harmony", dims = 1:50)
scRNA <- FindClusters(scRNA, resolution = 0.6)

Idents(scRNA) <- "seurat_clusters"
fimarker <- c(
  'CD3D',  
  "CD4","CD8A",
  "GPR183", 
  "FOXP3", 
  "CCR7","SELL",  
  "NKG7","GNLY", 
  "GZMB","GZMA","PRF1",
  "MKI67","STMN1", #Cycling
  "TRGV9","TRDV2",        
  "TRAV1-2"              
)
m <-DotPlot(scRNA,features = fimarker)+coord_flip()+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"), 
        axis.text.x=element_text(angle=60,hjust = 0.9,vjust=0.9))+
  scale_color_gradient(low="white",high="#8F4F17")+
  labs(x=NULL,y=NULL)
ggsave("scRNA_geneDotplot.pdf", m, width = 8, height = 8)

Idents(scRNA) <- "seurat_clusters"
celltype <- c("0"="TC", 
              "1"="TC", 
              "2"="mono", 
              "3"="TC",
              "4"="TC", 
              "5"="NK", 
              "6"="NK",
              "7"="BC",
              "8"="mono",
              "9"="BC", 
              "10"="TC", 
              "11"="TC", 
              "12"="TC",
              "13"="mono",
              "14"="mono", 
              "15"="TC", 
              "16"="DC",
              "18"="TC",
              "19"="cycling",
              "20"="DC",
              "21"="BC",
              "22"="progenitor",
              "23"="TC",
              "24"="mono",
              "25"="progenitor",
              "26"="BC"
)
scRNA <- RenameIdents(scRNA, celltype) 
scRNA$celltype <- scRNA@active.ident
scRNA$celltype_all <- factor(scRNA$celltype_all,levels=c("TC","NK","BC","mono","DC","progenitor","cycling"))

Idents(scRNA) <- "celltype"
q <- DimPlot(scRNA, reduction = "umap",label=F, raster=FALSE, cols = paletteDiscrete(values = unique(scRNA@meta.data$celltype)))+
  NoLegend()+
  labs(x = "UMAP1", y = "UMAP2", title = "scRNA") +
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
ggsave("UMAP_scRNA_celltype.pdf", q, width = 7, height = 7)

Idents(scRNA) <- "samples"
q <- DimPlot(scRNA, reduction = "umap",label=F, raster=FALSE, cols = paletteDiscrete(values = unique(scRNA@meta.data$samples)))+
  #NoLegend()+
  labs(x = "UMAP1", y = "UMAP2", title = "scRNA-seq") +
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
ggsave("UMAP_scRNA_samples.pdf", q, width = 7, height = 7)

#QC
pdf("scRNA_QC.pdf", width = 3, height = 6)
VlnPlot(scRNA, features = c("nFeature_RNA"),pt.size = 0,cols = "#191970")+ geom_boxplot(width=0.5,col = "black", fill = "white")
VlnPlot(scRNA, features = c("nCount_RNA"),pt.size = 0,cols = "#191970")+ geom_boxplot(width=0.5,col = "black", fill = "white")
VlnPlot(scRNA, features = c("percent.mt"),pt.size = 0,cols = "#191970")+ geom_boxplot(width=0.5,col = "black", fill = "white")
dev.off()

#
table(scRNA@meta.data$samples)
table(scRNA@meta.data$Clusters)
cM <- confusionMatrix(paste0(scRNA@meta.data$Clusters), paste0(scRNA@meta.data$samples))
cM <- data.frame(cM / Matrix::rowSums(cM))
range(cM) 
custom_order <- c("1", "2", "3", "4","5","6","7","8","9","10","11",
                  "12","13","14","15","16","17","18","19","20","21") 
cM <- cM[custom_order, ]

#library(ComplexHeatmap)
#library(circlize)
col_fun <- colorRamp2(
  c(0, 0.05, 0.1), 
  c("#8c510a", "white", "#01665e")
)
pdf("scRNA_prop_.pdf")
Heatmap(cM, name = "Fraction of cells",  
        cluster_rows = FALSE, 
        cluster_columns = FALSE
        )
dev.off()

#Cell Proportion
pplist = list()
sce_groups = c("TC","NK","BC","mono","DC","progenitor","cycling")
for(group_ in sce_groups){
  cellper_  = cellper %>% select(one_of(c('sample','group',group_)))
  colnames(cellper_) = c('sample','group','percent')
  cellper_$percent = as.numeric(cellper_$percent)
  cellper_ <- cellper_ %>% group_by(group) %>% mutate(upper =  quantile(percent, 0.75), 
                                                      lower = quantile(percent, 0.25),
                                                      mean = mean(percent),
                                                      median = median(percent))
  print(group_)
  print(cellper_$median)
  
  pp1 = ggplot(cellper_,aes(x=group,y=percent)) + 
    geom_jitter(width = 0.25,size = 3.5) +
    stat_summary(fun=mean, geom="point", color="grey40") +
    theme_cowplot() +  scale_x_discrete(limits = c("F_Healthy", "GO")) +
    theme(axis.text = element_text(size = 12),axis.title = element_text(size = 12),legend.text = element_text(size = 13),
          legend.title = element_text(size = 12),plot.title = element_text(size = 12,face = 'plain'),legend.position = 'none') + 
    labs(title = group_,y='Percentage') +
    geom_errorbar(aes(ymin = lower, ymax = upper),col = "grey40",width =  0.5)
  
  labely = max(cellper_$percent)
  compare_means(percent ~ group,  data = cellper_)
  my_comparisons <- list( c("F_Healthy", "GO") )
  pp1 = pp1 + stat_compare_means(comparisons = my_comparisons, label = "p.signif",size = 3.5,method = "wilxcon.test")
  pplist[[group_]] = pp1
}
library(cowplot)
p<-plot_grid(pplist[['TC']],pplist[['NK']],pplist[['BC']],
             pplist[['mono']],pplist[['DC']],pplist[['progenitor']],
             pplist[['cycling']])
ggsave("Cell Proportion.pdf",p, width = 18, height = 18, onefile=F)

write.csv(Cellratio,"Cellratio.csv")
Cellratio$cg <- paste(Cellratio$celltype, Cellratio$group, sep = "_")

p <- ggplot(Cellratio, aes(x = cg, y = Freq, fill = group, order = TRUE)) +
  geom_bar(stat = "summary", fun = "mean", position = position_dodge2(width = 0.9)) +  
  geom_errorbar(stat = "summary", fun.data = mean_se, position = position_dodge2(width = 0.9), width = 0.2) +  
  labs(title = "scRNA", x = NULL, y = "Fraction of cells") +  
  geom_point(aes(fill=sample),color = "black", size = 1) + 
  theme_minimal() +   
  theme(axis.text.x = element_text(angle = 60, hjust = 1),  

# save
ggsave("scRNA_Cell Proportion.pdf", p, width = 7, height = 5)


#cluster markers
Idents(scRNA) <- "celltype"
allmarker <-FindAllMarkers(scRNA, only.pos =T, min.pct = 0.5, logfc.threshold = 0.5)
markers_by_cluster <- split(allmarker, allmarker$cluster)
for (cluster_name in names(markers_by_cluster)) {
  cluster_data <- markers_by_cluster[[cluster_name]]
  file_name <- paste0(cluster_name, "_markers.csv")
  write.csv(cluster_data, file_name, row.names = FALSE)
}

Idents(scRNA) <- "celltype_all"
cellmarker <- FindAllMarkers(scRNA, only.pos = TRUE, min.pct = 0.05, logfc.threshold = 0.1)
top10 <- cellmarker %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
AverageExp <- AverageExpression(scRNA)
expr <- AverageExp$RNA
expr_normalize <- log2(expr+1)
top10_genelist <- top10$gene
expr_normalize_filted_1 <- expr_normalize[top10_genelist,]
write.csv(expr_normalize,"expression.csv")

annotation_col = data.frame(CellType = factor(rep(c("TC","NK","BC","mono","DC","progenitor","cycling"), 1)))
rownames(annotation_col) = c("TC","NK","BC","mono","DC","progenitor","cycling")
ann_colors = list(
  CellType = c(TC = "#1E1F57", NK = "#C8051D", BC = "#74117C", mono = "#1D7A33", DC = "#778BC7", progenitor = "#FDE20B", cycling = "#EE6722" )
)
p <- pheatmap(expr_normalize_filted_1,scale = "row", border = F,cluster_rows = FALSE,
         cluster_cols = FALSE,show_rownames = TRUE,  angle_col = 45, cellwidth = 25,
         annotation_col = annotation_col, annotation_legend = FALSE, annotation_colors = ann_colors, 
         color = colorRampPalette(colors = c("navy", "white", "firebrick3"))(100))
ggsave("scRNA_cellmarker_top10.pdf",p,width = 6, height = 9)

#DEG
cell_types <- rownames(table(scRNA$celltype))

for (cell_type in cell_types) {
  Idents(scRNA) <- "celltype"
  subset_data <- subset(scRNA, idents = c(cell_type))
  Idents(subset_data) <- "Diseasestate"
  subset_data <- subset(x = subset_data, downsample = 200)
  result <- FindMarkers(
    subset_data,
    ident.1 = "GO",
    ident.2 = "Healthy",
    group.by = "Diseasestate",
    logfc.threshold = 0.25,
    test.use = "MAST"
  )
  result <- result[result$p_val_adj < 0.05, ]
  file_name <- paste0(output_path, cell_type, "_DEG.csv")
  write.csv(result, file_name, row.names = TRUE)
}

list.files(output_path, pattern = "*.csv")

#
Idents(scRNA) <- "celltype"
TC <- subset(scRNA, idents = c("TC"))
NK <- subset(scRNA, idents = c("NK"))
BC <- subset(scRNA, idents = c("BC"))
mono <- subset(scRNA, idents = c("mono"))
DC <- subset(scRNA, idents = c("DC"))


