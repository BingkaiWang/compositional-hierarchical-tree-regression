rm(list=ls())
set.seed(123)
library(tidyverse)
library(xtable)
library(data.tree)
library(networkD3)
setwd("~/Dropbox/research/graphical-model/hierarchical-lasso/compositional-hierarchical-tree-regression/data-analysis/")
load("ADNI/MRICloud/200925/Data_vol_v1_Type1.RData")
thresholding <- function(x, threshold) {ifelse(abs(x) > threshold, x, 0)}
tuning_param_all <- readRDS("ADNI-CTASSO-BIC-allgroup.rds")[[1]]
alpha_allgroup <- readRDS("ADNI-CTASSO-BIC-allgroup.rds")[[2]]
alpha_allgroup <- round(alpha_allgroup - mean(alpha_allgroup),2)
beta_allgroup <- readRDS("ADNI-CTASSO-BIC-allgroup.rds")[[3]] %>% round(2)
my.roi <- readRDS("tree-structure.rds")

# tree visualization ==========
q <- ncol(my.roi)
my.roi$All_brain <- "Whole brain"
my.roi$pathString <- "Whole brain"
for(j in 1:nrow(my.roi)){
  string <- my.roi$pathString[j]
  for(k in (q-1):1){
    if(my.roi[j,k] != my.roi[j,k+1]){
      string <- paste(string, my.roi[j,k], sep = "/")
    }
  } 
  my.roi$pathString[j] <- gsub('.lvl[0-9]', '', string)
}
tree.structure <- as.Node(my.roi)
radialNetwork(ToListExplicit(tree.structure, unname = TRUE), fontSize = 9, opacity = 0.9, width = 800, height=800)


# alpha visualization ======== 
alpha_matrix <- data.frame(names = as.character(my.roi[,1]), alpha = alpha_allgroup[,1]) %>% filter(abs(alpha)>0)
alpha_matrix <- alpha_matrix[order(abs(alpha_matrix$alpha), decreasing = T),]
top10effect <- alpha_matrix[1:11,]
sum(abs(top10effect$alpha))/sum(abs(alpha_matrix$alpha))
write.table(data.frame(alpha = alpha_allgroup[,1], ROI.info$Index), "alpha-allgroup.txt", col.names = F)
alpha_matrix$names <- alpha_matrix$names %>% gsub('.lvl[0-9]', '', .) %>% gsub('_', '-', .)
print(xtable(cbind(alpha_matrix[1:17,], alpha_matrix[18:34,], alpha_matrix[c(35:50, 50),])), include.rownames=FALSE)

# beta visualization=========
beta_alpha_top10effect <- beta_allgroup[my.roi[my.roi$level5 %in% top10effect$names, 1:5] %>% unlist %>% unique,]
beta_alpha_top10effect <- thresholding(beta_alpha_top10effect, 0)
c_scale <- colorRamp(c('blue', 'white', 'red'))
cbind(names(beta_alpha_top10effect), round(c_scale((beta_alpha_top10effect/358 + 1)/2)))

betaCSF <- data.frame(beta=ifelse(my.roi$level1 == "CSF.lvl1" & my.roi$level5 %in% top10effect$names, alpha_allgroup, 0), ROI.info$Index)
write.table(betaCSF, "beta1.txt", col.names = F)
cbind(betaCSF[abs(betaCSF$beta) > 0, 1], round(c_scale((betaCSF[abs(betaCSF$beta) > 0, 1]/540 + 1)/2)))
betaTL <- data.frame(beta=ifelse(my.roi$level1 == "Telencephalon_L.lvl1" & my.roi$level5 %in% top10effect$names, alpha_allgroup, 0), ROI.info$Index)
write.table(betaTL, "beta2.txt", col.names = F)
cbind(betaTL[abs(betaTL$beta) > 0, 1], round(c_scale((betaTL[abs(betaTL$beta) > 0, 1]/540 + 1)/2)))

beta_matrix <- data.frame(names = names(beta_allgroup[abs(round(beta_allgroup))>0,]),
                          beta = beta_allgroup[abs(round(beta_allgroup))>0,]) %>% arrange(1/abs(beta))
print(xtable(cbind(beta_matrix[1:32,], beta_matrix[33:64,], beta_matrix[c(65:95, 95),])), include.rownames=FALSE)
