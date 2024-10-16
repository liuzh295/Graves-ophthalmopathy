library(Seurat)
library(scRepertoire)
library(dplyr)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(viridis)
library(DoubletFinder)
library(ArchR)
library(ggpubr)
library(readxl)
library(tidyverse)
library(tidyr)
library(ggplot2)
library(dplyr)
library(scales)

contig_list <- list(H_02,H_03,H_04,H_05,H_06,H_07,H_08,H_09,H_10,H_11,
                    GO_01,GO_02,GO_03,GO_04,GO_05,GO_06,GO_07,GO_08,GO_09,GO_10,GO_11,
                    GO_12,GO_13,GO_14,GO_15,GO_16,GO_17,GO_18,GO_19,GO_20,GO_21)
TCR <-combineTCR(
  contig_list,
  samples = c("H_02","H_03","H_04","H_05","H_06","H_07","H_08","H_09","H_10","H_11",
              "GO_01","GO_02","GO_03","GO_04","GO_05","GO_06","GO_07","GO_08","GO_09","GO_10",
              "GO_11","GO_12","GO_13","GO_14","GO_15","GO_16","GO_17","GO_18","GO_19","GO_20","GO_21"), 
  ID = c("Healthy","Healthy","Healthy","Healthy","Healthy","Healthy","Healthy","Healthy","Healthy","Healthy",
         "GO","GO","GO","GO","GO","GO","GO","GO","GO","GO","GO","GO","GO","GO","GO","GO","GO","GO","GO","GO","GO"),
  #cells = c("T-AB"),
  removeNA = FALSE,
  filterMulti = FALSE,
  removeMulti = FALSE
)

scTCR <- combineExpression(TCR, obj, 
                          cloneCall="gene", 
                          #group.by = "samples", 
                          proportion = FALSE, 
                          cloneTypes=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500), 
                          filterNA = FALSE)


#diversity
clonalDiversity(TCR, 
                cloneCall = "gene", 
                group.by = "sample",
                x.axis = "ID",
                n.boots = 100)

#uniqueclone
quantContig(TCR, cloneCall="gene+nt", scale = F)

#clonetype
ggplot(uniqueGO, aes(x = celltypeTC, fill = cloneType)) +
  geom_bar(position = "fill") + 
  scale_fill_manual(values = rep(c("#FF4B20", "#FFB433", "#C6FDEC", "#7AC5FF", "#0347A6"),13)) +  
  scale_y_continuous(labels = percent) +    
  guides(fill=guide_legend(title = "cloneType")) +
  labs(title = "TCR",
       x = " ",
       y = "percent of GO specific clonetypes") +
  theme_minimal()+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"), 
        axis.text.x=element_text(angle=60,hjust = 0.9,vjust=0.9))
ggsave("TCR_GOspecificclonetype_percent.pdf", width = 6, height = 6)


