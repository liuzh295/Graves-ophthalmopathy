library(ArchR)
library(Seurat)
library(tidyverse)
addArchRGenome("hg38")
addArchRThreads(threads = 30)

Idents(scRNA) <- "celltype"
addArchRThreads(threads = 1)
Integrate <- addGeneIntegrationMatrix(
  ArchRProj = scATAC, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  scRNA = scRNA,
  addToArrow = FALSE,
  force= TRUE,
  #groupList = groupList,
  groupRNA = "celltype",
  nameCell = "predictedCell",
  nameGroup = "predictedGroup",
  nameScore = "predictedScore"
)
cM <- as.matrix(confusionMatrix(Integrate$Clusters, Integrate$predictedGroup))
preClust <- colnames(cM)[apply(cM, 1 , which.max)]
cbind(preClust, rownames(cM)) #Assignments
unique(unique(Integrate$predictedGroup))
pal <- paletteDiscrete(values = scRNA@meta.data$celltype)
p1 <- plotEmbedding(
  Integrate, 
  embedding = "UMAPHarmony",
  colorBy = "cellColData", 
  name = "predictedGroup", 
  pal = pal
)
plotPDF(p1, name = "UMAP-scATAC+RNA-Integration_predictedGroup_.pdf", 
        ArchRProj = Integrate, addDOC = FALSE, width = 5, height = 5)

#
library(BSgenome.Hsapiens.UCSC.hg38)
Integrate <- addGroupCoverages(ArchRProj = Integrate, groupBy = "predictedGroup")
#####peakcalling
pathToMacs2 <- findMacs2()
Integrate <- addReproduciblePeakSet(
  ArchRProj = Integrate, 
  groupBy = "predictedGroup", 
  pathToMacs2 = pathToMacs2
)
peakset <- getPeakSet(Integrate)
Integrate <- addPeakMatrix(Integrate)
Integrate <- addMotifAnnotations(ArchRProj = Integrate, motifSet = "cisbp", name = "Motif")
names(Integrate@peakAnnotation)
Integrate <- addBgdPeaks(Integrate)
Integrate <- addDeviationsMatrix(
  ArchRProj = Integrate, 
  peakAnnotation = "Motif",
  force = TRUE
)

#
#########Peak2GeneLinkages
Integrate <- addPeak2GeneLinks(
  ArchRProj = Integrate,
  reducedDims = "IterativeLSI"
)
p2g <- getPeak2GeneLinks(
  ArchRProj = Integrate,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)
p1 <- plotPeak2GeneHeatmap(ArchRProj = Integrate, groupBy = "predictedGroup",k = 27)
plotPDF(p1, name = "UMAP-Integration_peak2link.pdf", 
        ArchRProj = Integrate, addDOC = FALSE, width = 7, height = 10)

#
#########TF
seGroupMotif <- getGroupSE(ArchRProj = Integrate, useMatrix = "MotifMatrix", groupBy = "predictedGroup")
seZ <- seGroupMotif[rowData(seGroupMotif)$seqnames=="z",]
rowData(seZ)$maxDelta <- lapply(seq_len(ncol(seZ)), function(x){
  rowMaxs(assay(seZ) - assay(seZ)[,x])
}) %>% Reduce("cbind", .) %>% rowMaxs

corGSM_MM <- correlateMatrices(
  ArchRProj = Integrate,
  useMatrix1 = "GeneScoreMatrix",
  useMatrix2 = "MotifMatrix",
  reducedDims = "IterativeLSI"
)

corGSM_MM$maxDelta <- rowData(seZ)[match(corGSM_MM$MotifMatrix_name, rowData(seZ)$name), "maxDelta"]
corGSM_MM <- corGSM_MM[order(abs(corGSM_MM$cor), decreasing = TRUE), ]
corGSM_MM <- corGSM_MM[which(!duplicated(gsub("\\-.*","",corGSM_MM[,"MotifMatrix_name"]))), ]
corGSM_MM$TFRegulator <- "NO"
corGSM_MM$TFRegulator[which(corGSM_MM$cor > 0.5 & corGSM_MM$padj < 0.01 & corGSM_MM$maxDelta > quantile(corGSM_MM$maxDelta, 0.75))] <- "YES"
sort(corGSM_MM[corGSM_MM$TFRegulator=="YES",1])
p <- ggplot(data.frame(corGSM_MM), aes(cor, maxDelta, color = TFRegulator)) +
  geom_point() + 
  theme_ArchR() +
  geom_vline(xintercept = 0, lty = "dashed") + 
  scale_color_manual(values = c("NO"="darkgrey", "YES"="firebrick3")) +
  xlab("Correlation To Gene Score") +
  ylab("Max TF Motif Delta") +
  scale_y_continuous(
    expand = c(0,0), 
    limits = c(0, max(corGSM_MM$maxDelta)*1.05)
  )
ggsave("Genescore_TF.pdf",p)
data <- data.frame(corGSM_MM)
write.csv(data, "Genescore_TF.csv")

corGIM_MM <- correlateMatrices(
  ArchRProj = Integrate,
  useMatrix1 = "GeneIntegrationMatrix",
  useMatrix2 = "MotifMatrix",
  reducedDims = "IterativeLSI"
)
corGSM_MM$maxDelta <- rowData(seZ)[match(corGSM_MM$MotifMatrix_name, rowData(seZ)$name), "maxDelta"]
corGSM_MM <- corGSM_MM[order(abs(corGSM_MM$cor), decreasing = TRUE), ]
corGSM_MM <- corGSM_MM[which(!duplicated(gsub("\\-.*","",corGSM_MM[,"MotifMatrix_name"]))), ]
corGSM_MM$TFRegulator <- "NO"
corGSM_MM$TFRegulator[which(corGSM_MM$cor > 0.5 & corGSM_MM$padj < 0.01 & corGSM_MM$maxDelta > quantile(corGSM_MM$maxDelta, 0.75))] <- "YES"
sort(corGSM_MM[corGSM_MM$TFRegulator=="YES",1])

p <- ggplot(data.frame(corGSM_MM), aes(cor, maxDelta, color = TFRegulator)) +
  geom_point() + 
  theme_ArchR() +
  geom_vline(xintercept = 0, lty = "dashed") + 
  scale_color_manual(values = c("NO"="darkgrey", "YES"="firebrick3")) +
  xlab("Correlation To Gene Score") +
  ylab("Max TF Motif Delta") +
  scale_y_continuous(
    expand = c(0,0), 
    limits = c(0, max(corGSM_MM$maxDelta)*1.05)
  )
ggsave("Geneinte_TF_.pdf",p)
data <- data.frame(corGSM_MM)
write.csv(data, "Geneinte_TF_.csv")
