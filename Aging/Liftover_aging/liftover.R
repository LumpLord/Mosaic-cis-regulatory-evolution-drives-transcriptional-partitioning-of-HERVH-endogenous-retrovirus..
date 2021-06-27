---
title: "ltr7LO"
author: "me"
date: "2/12/2021"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```



```{r cars}
library(ggplot2)
library(RColorBrewer)


setwd("C:/Users/thoma/Documents/LTR_recombination/7up/figures/3")

LO=read.csv("LO_gross.csv", row.names=1)
#LO=data.frame("liftover.csv", row.names=1)
LO
cols <- c(up1='darkred', up2='chocolate', u1='lightpink' , u2='green' , d1='magenta' , d2='purple1' , bc='cyan' , o='royalblue2' , b='gold' , c='lightskyblue2' , y='seagreen')

p<-ggplot(LO, aes(x=MYA, y=LO, group=Subfamily)) +
  geom_line(aes(color=Subfamily))+
  geom_point(aes(color=Subfamily))#+
  #scale_color_manual(cols)
  #scale_color_manual('darkred','chocolate', 'lightpink' , 'green' , 'magenta' ,'purple1' , 'cyan' , 'royalblue2' ,  'seagreen' , 'gold' , 'lightskyblue2') 

p+scale_color_manual(values=cols)+theme_light()
ggsave("7LO_gross.pdf", width = 15, height = 10)

```