---
title: "7subs_heat"
author: "me"
date: "8/28/2020"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## R Markdown
```{r}
library(pheatmap)
library(RColorBrewer)
library(viridis)

setwd("C:/Users/thoma/Documents/LTR_recombination/7up/figures/2")


heat=read.csv("heat_7uptranscribed2.csv", row.names=1)
#heat=as.data.frame(t(scale(t("heatmap_props_full.csv"), scale=F)))
rows <- rownames(heat)
heat
#heated <- as.data.frame(sapply(heat, as.character))
#heated <- as.data.frame(t(scale(t(heat), scale=F)))
#rownames(heated) <- rows
#heated

#heatmap.colours <- c("white","grey","seagreen3","darkgreen",
                    #"green","brown","tan","red","orange",
                    #"pink","magenta","purple","blue","skyblue3",
                    #"blue","skyblue2")
heat.colours <- c(viridis(100))
#heat.colours <- colorRampPalette(brewer.pal(9,"Blues"))(100)

#names(heat.colours) <- 1

# The mtcars dataset:
#heatboi <- as.matrix(heat)

# Default Heatmap
heatmap(heatboi, col = topo.colors(50))


#gheatmap(data=heat, offset = 0.01, color=NULL, 
         #colnames_position="top", 
         #colnames_angle=90, colnames_offset_y = 1, 
         #hjust=0, font.size=2) 

#ggplot(heat, aes(X, Y, fill= Z)) + 
  #geom_tile() +
  #scale_fill_grad

#ggsave("Heat_full.pdf", width = 10, height = 10)

p <- pheatmap(
    mat               = heat,
  main              = "7up1/2 Highly Transcribed",
  color             = heat.colours,
  border_color      = NA,
  show_colnames     = TRUE,
  show_rownames     = TRUE,
#  annotation_col    = col_annot,
  drop_levels       = TRUE,
#   cluster_col    = FALSE,
#  annotation_names_row = F,
  fontsize_row          = 10,
  fontsize_col = 8,
  filename = "Heat_7uptrans_limited.pdf"
)
  
```

