library(ArchR)
library(Seurat)
library(tidyverse)
addArchRGenome("hg38")
addArchRThreads(threads = 20)

scATAC <- addIterativeLSI(
  ArchRProj = scATAC,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  iterations = 2, 
  clusterParams = list(
    resolution = c(0.2), 
    sampleCells = 10000, 
    n.start = 10
  ), 
  varFeatures = 15000, 
  dimsToUse = 1:30,
  force = TRUE
)

scATAC <- addHarmony(
  ArchRProj = scATAC,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "Sample", 
  force = TRUE
)

scATAC <- addClusters(
  input = scATAC,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.6, #resolution
  sampleCells = 50000, 
  force = TRUE
)

scATAC <- addUMAP(
  ArchRProj = scATAC, 
  reducedDims = "Harmony", 
  name = "UMAPHarmony", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine", 
  force = TRUE
)
p1 <- plotEmbedding(ArchRProj = scATAC, colorBy = "cellColData", name = "Sample", embedding = "UMAPHarmony")
p2 <- plotEmbedding(ArchRProj = scATAC, colorBy = "cellColData", name = "Clusters", embedding = "UMAPHarmony")
plotPDF(p1, name = "UMAP-Harmony-Sample.pdf", ArchRProj = scATAC, addDOC = FALSE, width =5, height = 5)
plotPDF(p2, name = "UMAP-Harmony-Clusters.pdf", ArchRProj = scATAC, addDOC = FALSE, width =5, height = 5)

markersGS <- getMarkerFeatures(
  ArchRProj = scATAC, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
markerGenes  <- c(
  "markers"
)

plotPDF(plotList = p, 
        name = "Plot-UMAP-Marker-Genes-W-Imputation.pdf", 
        ArchRProj = scATAC, 
        addDOC = FALSE, width = 5, height = 5)

p <- plotBrowserTrack(
  ArchRProj = scATAC, 
  groupBy = "Clusters", 
  geneSymbol = markerGenes, 
  upstream = 50000,
  downstream = 50000
)
plotPDF(plotList = p, 
        name = "Plot-Tracks-Marker-Genes_magic.pdf", 
        ArchRProj = scATAC, 
        addDOC = FALSE, width = 5, height = 5)

cM <- confusionMatrix(scATAC$Clusters, scATAC$predictedGroup)
labelOld <- rownames(cM)
labelOld
remapClust <- c(
  "0"="xxx", 
  "1"="xxx", 
  "2"="xxx", 
  "3"="xxx",
  "4"="xxx", 
  "5"="xxx", 
  "6"="xxx"
)
scATAC$Clusters2 <- mapLabels(scATAC$Clusters, newLabels = remapClust, oldLabels = labelOld)
p1 <- plotEmbedding(scATAC, colorBy = "cellColData", name = "Clusters2")
p1
plotPDF(p1, name = "Plot-UMAP-Remap-Clusters.pdf", ArchRProj = projHeme2, addDOC = FALSE, width = 5, height = 5)

#Cell Proportion
table(scATAC$Sample)
table(scATAC$celltype_all, scATAC$Sample)
Cellratio <- prop.table(table(scATAC$celltype_all, scATAC$Sample), margin = 2)
Cellratio <- data.frame(Cellratio)
write.csv(Cellratio,"Cellratio.csv")
Cellratio <- read.csv("Cellratio_.csv")
Cellratio$cg <- paste(Cellratio$celltype, Cellratio$group, sep = "_")
p <- ggplot(Cellratio, aes(x = cg, y = Freq, fill = group, order = TRUE)) +
  geom_bar(stat = "summary", fun = "mean", position = position_dodge2(width = 0.9)) + 
  geom_errorbar(stat = "summary", fun.data = mean_se, position = position_dodge2(width = 0.9), width = 0.2) +  
  labs(title = "scATAC", x = NULL, y = "Fraction of cells") +  
  geom_point(aes(fill=sample),color = "black", size = 1) +  
  theme_minimal() +   
  theme(axis.text.x = element_text(angle = 60, hjust = 1),  
        panel.grid = element_blank())
ggsave("scATAC_Cell Proportion.pdf", p, width = 7, height = 5)

#cluster markers
scATAC <- readRDS("scATAC_scATAC_final.rds")
table(scATAC@cellColData$celltype)
markersPeaks <- getMarkerFeatures(
  ArchRProj = scATAC, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "celltype",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon")
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
lapply(names(markerList), function(cell_name) {
  cell_data <- markerList[[cell_name]]
  file_name <- paste0(gsub(" ", "_", cell_name), "_markers.csv")
  write.csv(cell_data, file_name, row.names = FALSE)
})

#DAR
table(scATAC@cellColData$celltype)
cell_types <- table(scATAC$celltype)
selected_cells <- c()
for (cell_type in names(cell_types)) {
  idxSample <- BiocGenerics::which(scATAC$celltypescATAC == cell_type)
  if (length(idxSample) >= 1000) {
    sampled_cells <- sample(scATAC$cellNames[idxSample], size = 1000, replace = FALSE)
  } else {
    sampled_cells <- scATAC$cellNames[idxSample]
  }
  selected_cells <- c(selected_cells, sampled_cells)
}
scATAC <- scATAC[selected_cells, ]
table(scATAC@cellColData$DiseaseStates)
table(scATAC$celltypescATAC)
for (cell_type in cell_types) {
  idxSample <- BiocGenerics::which(scATAC$celltypescATAC %in% c(cell_type))
  cellsSample <- scATAC$cellNames[idxSample]
  scATAC <- scATAC[cellsSample, ]
  cat("Processing:", cell_type, "\n")
  table(scATAC$celltypescATAC)
  markersPeaks <- getMarkerFeatures(
    ArchRProj = scATAC, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "DiseaseStates",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon",
    useGroups = "GO",
    bgdGroups = "F-H"
  )
  markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 0.5")
  data <- markerList$"GO"
  output_file <- paste0(cell_type, "_DAR.csv")
  write.csv(data, output_file)
}


