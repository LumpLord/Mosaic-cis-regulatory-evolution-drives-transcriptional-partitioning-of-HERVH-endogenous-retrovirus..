---
title: "TBL_7up"
author: "me"
date: "8/4/2020"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_knit$set(root.dir = "C:/Users/thoma/Documents/LTR_recombination/7up/newtrees/LTR7only_final/adj/group_updates")
```


```{r}
# Libraries
library(tidyverse)
library(hrbrthemes)
library(viridis)
library(ggplot2)
library(dplyr)
```

```{r}
#set wd
setwd("C:/Users/thoma/Documents/LTR_recombination/7up/newtrees/LTR7only_final/adj/group_updates")
```


```{r}
# load a dataset
dat1 = read.csv("7sub_tbl.csv")
dat1$subfamily <- factor(dat1$subfamily, levels = c("7o", "7bc", "7d1", "7d2", "7u1", "7u2", "7up2", "7up1"))
dat1
```


```{r}

# Plot
dat1_graph =  ggplot(data = dat1, aes(x=subfamily, y=terminal.branch.length, fill=subfamily))


dat1_graph +  geom_boxplot() +
    scale_fill_viridis(discrete = TRUE, alpha=0.6) +
    geom_jitter(color="black", size=0.4, alpha=0.9) +
    theme_ipsum() +
    theme(
      legend.position="none",
      plot.title = element_text(size=11)
    ) +
    ggtitle("A boxplot with jitter") +
    xlab("") + ggsave("7sub_tbl.pdf", device="pdf", scale=1, dpi=300)

dat1_graph
```
