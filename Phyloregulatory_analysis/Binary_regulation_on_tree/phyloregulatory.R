---
title: "7_circ"
author: "me"
date: "8/26/2020"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## R Markdown

```{r}
setwd("C:/Users/thoma/Documents/LTR_recombination/7up/figures/circTree/bigtree")
#library(phyloseq)
library(ggjoy)
library(dplyr)
library(ggtree)
library(RColorBrewer)

info <- read.csv("bigtreematrix.csv")
tree <- read.tree("LTR7only_final_groups.nwk")
cols <- c(up1='darkred', up2='chocolate', u1='gold' , u2='green' , d1='magenta' , d2='purple1' , bc='cyan' , o='royalblue2')

p <- ggtree(tree, layout='circular') %<+% info + 
  geom_tippoint(aes(color=subfamily),size=1.0) + 
  scale_color_manual(values=cols) #+ 
  #geom_tiplab2(aes(label=Transcribed), align=T, linetype=NA, 
              #size=2, offset=4, hjust=0.5)# +
  #geom_tiplab2(aes(Top 7th transcribed), align=T, linetype=NA, 
              #size=2, offset=8, hjust=0.5)


heatmapData=read.csv("bigtreematrix_heat.csv", row.names=1)
rn <- rownames(heatmapData)
heatmapData <- as.data.frame(sapply(heatmapData, as.character))
rownames(heatmapData) <- rn

#heatmap.colours <- c("white","grey","seagreen3","darkgreen",
                    #"green","brown","tan","red","orange",
                    #"pink","magenta","purple","blue","skyblue3",
                    #"blue","skyblue2")
heatmap.colours <- c("white", "set3")
names(heatmap.colours) <- 0:1

q <- open_tree(p, angle=15)
q

gheatmap(q, heatmapData, offset = 0.01, color=NULL, 
         colnames_position="top", 
         colnames_angle=90, colnames_offset_y = 1, 
         hjust=0, font.size=4) 

#gheatmap(p, heatmapData, offset = 10, color=NULL, 
        # colnames_position="top", 
         #colnames_angle=90, colnames_offset_y = 1, 
         #hjust=0, font.size=2) #+
  #scale_fill_manual(values=heatmap.colours, breaks=0:15)

#fasta <- system.file("7up_newphylo_names.fas", package="ggtree")
#msaplot(ggtree(beast_tree), fasta) 



ggsave("7all_regTree_reorder_open.pdf", width = 20, height = 20)
```

