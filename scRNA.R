library(dplyr)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(viridis)
library(ArchR)
library(harmony)
library(ggpubr)

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
fimarker <- c("markers"
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
celltype <- c("0"="xxx", 
              "1"="xxx", 
              "2"="xxx", 
              "3"="xxx",
              "4"="xxx", 
              "5"="xxx", 
              "6"="xxx"
)
scRNA <- RenameIdents(scRNA, celltype) 
scRNA$celltype <- scRNA@active.ident

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

#Cell Proportion
pplist = list()
sce_groups = c("xxx", "xxx", "xxx", "xxx")
library(ggplot2)
library(dplyr)
library(ggpubr)
library(cowplot)
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
  geom_bar(stat = "summary", fun = "mean", position = position_dodge2(width = 0.9)) +  # 设置position为position_dodge2并调整width
  geom_errorbar(stat = "summary", fun.data = mean_se, position = position_dodge2(width = 0.9), width = 0.2) +  # 添加误差线
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

#DEG
cell_types <- rownames(table(seRNA$celltypeTC))

for (cell_type in cell_types) {
  Idents(seRNA) <- "celltype"
  subset_data <- subset(seRNA, idents = c(cell_type))
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




