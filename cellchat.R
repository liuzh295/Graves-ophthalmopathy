library(CellChat)
library(patchwork)
library(ggalluvial)
library(igraph)
library(dplyr)
options(stringsAsFactors = FALSE) 
options(futrue.globlas.Maxsize=2*1024**3)

colnames(TAO_scRNA@meta.data)
table(TAO_scRNA$Diseasestate)
Idents(TAO_scRNA) <- "Diseasestate"
sco.HD <- subset(TAO_scRNA, idents = c("F_Healthy"))
sco.GO <- subset(TAO_scRNA, idents = c("GO"))
###
cco.HD <- createCellChat(sco.HD@assays$RNA@data, meta = sco.HD@meta.data, group.by = "celltypeN")
cco.GO <- createCellChat(sco.GO@assays$RNA@data, meta = sco.GO@meta.data, group.by = "celltypeN")
save(cco.HD,cco.GO, file="cco.rds")
#
future::plan("multisession", workers = 4)
options(future.globals.maxSize= 8000*1024^2)
cellchat <- cco.HD
cellchat@DB <- CellChatDB.human
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, raw.use = T, population.size = T)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)#
cellchat <- netAnalysis_computeCentrality(cellchat,slot.name = "netP")
cco.HD <- cellchat
saveRDS(cco.HD,"cco.HD.rds")

future::plan("multisession", workers = 4)
options(future.globals.maxSize= 8000*1024^2)
cellchat <- cco.GO
cellchat@DB <- CellChatDB.human
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, raw.use = T, population.size = T)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)#
cellchat <- netAnalysis_computeCentrality(cellchat,slot.name = "netP")
cco.GO <- cellchat
saveRDS(cco.GO,"cco.GO.rds")

#
cco.list <- list(hd=cco.HD, go= cco.GO)
cellchat <- mergeCellChat(cco.list, add.names = names(cco.list), cell.prefix= TRUE)

#
g1 <- compareInteractions(cellchat,show.legend = F, group = c(1,2),measure = "count")
g2 <- compareInteractions(cellchat,show.legend = F, group = c(1,2),measure = "weight")
p <- g1+g2
ggsave("communicationcountandweight.pdf", p, width= 6, height = 4)


cco.list2 <- list(hd=cco.HD, go= cco.GO)
cellchat2 <- mergeCellChat(cco.list2, add.names = names(cco.list2), cell.prefix= TRUE)
pdf("diffInteractionGOVSH_heatmap.pdf")
netVisual_heatmap(cellchat2)
netVisual_heatmap(cellchat2, measure = "weight")
dev.off()
pdf("diffInteractionGOVSH.pdf")
netVisual_diffInteraction(cellchat2,weight.scale=T)
netVisual_diffInteraction(cellchat2,weight.scale=T,measure = "weight")
dev.off()

weight.max <- getMaxWeight(cco.list2, attribute= c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
pdf("weight.max.pdf")
for (i in 1:length(cco.list2)) {
  netVisual_circle(cco.list2[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(cco.list2)[i]))
}
dev.off()

gg1 <- rankNet(cellchat2, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat2, mode = "comparison", stacked = F, do.stat = TRUE)
p<-gg1 + gg2
ggsave("h_go.pdf", p, width= 10, height = 5)

library(ComplexHeatmap)
i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(cco.list2[[1]]@netP$pathways, cco.list2[[2]]@netP$pathways)
pdf("HCGO_all.pdf", width = 12, height = 7)
ht1 = netAnalysis_signalingRole_heatmap(cco.list[[1]], pattern = "all", signaling = pathway.union, title = names(cco.list)[1], width = 11, height = 12, color.heatmap = "Reds")
ht2 = netAnalysis_signalingRole_heatmap(cco.list[[2]], pattern = "all", signaling = pathway.union, title = names(cco.list)[2], width = 11, height = 12, color.heatmap = "Reds")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()
pdf("HCGO_outgoing.pdf", width = 12, height = 7)
ht1 = netAnalysis_signalingRole_heatmap(cco.list[[1]], pattern = "outgoing", signaling = pathway.union, title = names(cco.list)[1], width = 11, height = 12,color.heatmap = "YlOrRd")
ht2 = netAnalysis_signalingRole_heatmap(cco.list[[2]], pattern = "outgoing", signaling = pathway.union, title = names(cco.list)[2], width = 11, height = 12,color.heatmap = "YlOrRd")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()
pdf("HCGO_incoming.pdf", width = 12, height = 7)
ht1 = netAnalysis_signalingRole_heatmap(cco.list[[1]], pattern = "incoming", signaling = pathway.union, title = names(cco.list)[1], width = 11, height = 12,color.heatmap = "RdPu")
ht2 = netAnalysis_signalingRole_heatmap(cco.list[[2]], pattern = "incoming", signaling = pathway.union, title = names(cco.list)[2], width = 11, height = 12,color.heatmap = "RdPu")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

#CD8 Te
levels(cellchat@idents$joint)
pdf("tissueCD8Te.pdf", width = 12, height = 10) 
netVisual_bubble(cellchat, sources.use = 7, targets.use = c((1:6),(8:23)),  comparison = c(1, 2), angle.x = 45)
dev.off()

pdf("tissueCD8Te_.pdf", width = 12, height = 10) 
netVisual_bubble(cellchat, sources.use = c((1:6),(8:23)), targets.use = 7,  comparison = c(1, 2), angle.x = 45)
dev.off()
