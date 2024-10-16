library(Seurat)
library(patchwork)
library(clusterProfiler)
library(org.Hs.eg.db) 
library(ggplot2)
library(Seurat)
library(patchwork)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(forcats)

#DEG-GO
sample_names <- c("xxx","xxx","xxx")
ego_list <- list()
fraction_to_decimal <- function(fraction) {
  numerator <- as.numeric(gsub(".*?([0-9]+)/.*", "\\1", fraction))
  denominator <- as.numeric(gsub(".*?/([0-9]+)", "\\1", fraction))
  return(numerator / denominator)
}
for (sample in sample_names) {
  DEG <- read.csv(paste0(sample, ".csv"))
  ego <- enrichGO(gene = DEG$X,
                  OrgDb = "org.Hs.eg.db",
                  keyType = "SYMBOL",
                  ont = "ALL",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 1,
                  qvalueCutoff = 1)
  
  ego_df <- data.frame(ego)
  ego_df$Group <- sample
  ego_df$GeneRatio_ <- sapply(ego_df$GeneRatio, fraction_to_decimal)
  ego_df <- ego_df[order(ego_df$p.adjust, decreasing = FALSE), ]
  ego_df <- ego_df[1:5, ]
  ego_list[[sample]] <- ego_df
}
all <- do.call(rbind, ego_list)
all$Description <- as.factor(all$Description)
all$Description <- fct_inorder(all$Description)
all$value <- -log10(all$p.adjust)
ggplot(all, aes(Group, Description)) +
  geom_point(aes(color = value, size = GeneRatio_)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 0.9, vjust = 0.9)) +
  scale_color_gradient(low = '#6699CC', high = '#CC3333') +
  labs(x = NULL, y = NULL) +
  guides(size = guide_legend(order = 1)) +
  xlim(sample_names)
ggsave("xxx_DEGGO.pdf", width = 8, height = 8)

#DAR
sample_names <- c("xxx","xxx","xxx")
ego_list <- list()
fraction_to_decimal <- function(fraction) {
  numerator <- as.numeric(gsub(".*?([0-9]+)/.*", "\\1", fraction))
  denominator <- as.numeric(gsub(".*?/([0-9]+)", "\\1", fraction))
  return(numerator / denominator)
}
for (sample in sample_names) {
  DEG <- read.csv(paste0(sample, "_DAR.csv"))
  ego <- enrichGO(gene = DEG$name,
                  OrgDb = "org.Hs.eg.db",
                  keyType = "SYMBOL",
                  ont = "ALL",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 1,
                  qvalueCutoff = 1)
  
  ego_df <- data.frame(ego)
  ego_df$Group <- sample
  ego_df$GeneRatio_ <- sapply(ego_df$GeneRatio, fraction_to_decimal)
  ego_df <- ego_df[order(ego_df$p.adjust, decreasing = FALSE), ]
  ego_df <- ego_df[1:5, ]
  ego_list[[sample]] <- ego_df
}
all <- do.call(rbind, ego_list)
all$Description <- as.factor(all$Description)
all$Description <- fct_inorder(all$Description)
all$value <- -log10(all$p.adjust)
ggplot(all, aes(Group, Description)) +
  geom_point(aes(color = value, size = GeneRatio_)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 0.9, vjust = 0.9)) +
  scale_color_gradient(low='#1D7A33',high='#FDE20B') +
  labs(x = NULL, y = NULL) +
  guides(size = guide_legend(order = 1)) +
  xlim(sample_names)
ggsave("xxx_DARGO.pdf", width = 8, height = 8)


