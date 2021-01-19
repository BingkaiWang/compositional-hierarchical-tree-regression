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

# tree visualization
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


# alpha visualization
alpha_matrix <- data.frame(names = as.character(my.roi[,1]), alpha = alpha_allgroup[,1]) %>% filter(abs(alpha)>0)
alpha_matrix <- alpha_matrix[order(abs(alpha_matrix$alpha), decreasing = T),]
top10effect <- alpha_matrix[1:11,]
sum(abs(top10effect$alpha))/sum(abs(alpha_matrix$alpha))
write.table(data.frame(alpha = alpha_allgroup[,1], ROI.info$Index), "alpha-allgroup.txt", col.names = F)
alpha_matrix$names <- alpha_matrix$names %>% gsub('.lvl[0-9]', '', .) %>% gsub('_', '-', .)
print(xtable(cbind(alpha_matrix[1:17,], alpha_matrix[18:34,], alpha_matrix[c(35:50, 50),])), include.rownames=FALSE)

# beta visualization 
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


# beta_matrix <- matrix(NA, nrow = nrow(my.roi), ncol = ncol(my.roi)-1)
# for(j in 1: ncol(beta_matrix)){
#   for(i in 1: nrow(beta_matrix)){
#     beta_matrix[i,j] <- round(beta_allgroup[my.roi[i,j], 1], 0)
#   }
# }
# 
# beta_matrix[which(my.roi$level1 == "CSF.lvl1"),]

# alpha_matrix$names <- as.character(alpha_matrix$names)
# cbind(alpha_matrix[,2], my.roi[my.roi$level5 %in% alpha_matrix[,1], 1:4])
# alpha_matrix$names[alpha_matrix$names %in% c("STG_R.lvl5", "STG_R_pole.lvl5")] <- "STG_R.lvl4"
# alpha_matrix$names[alpha_matrix$names %in% c("STG_L.lvl5", "STG_L_pole.lvl5")] <- "STG_L.lvl4"
# alpha_matrix$names[alpha_matrix$names %in% c("LV_atrium_L.lvl5", "LV_Occipital_L.lvl5")] <- "PosteriorLateralVentricle_L.lvl4"
# alpha_matrix$names[alpha_matrix$names %in% c("LV_atrium_R.lvl5", "LV_Occipital_R.lvl5")] <- "PosteriorLateralVentricle_R.lvl4"
# alpha_matrix$names[alpha_matrix$names %in% c("SylFrontSul_L.lvl4", "SylTempSul_L.lvl4", " SylParieSul_L.lvl4")] <- "SylvianFissureExt_L.lvl3"
# alpha_matrix$names[alpha_matrix$names %in% c("SylFrontSul_R.lvl4", "SylTempSul_R.lvl4", " SylParieSul_R.lvl4")] <- "SylvianFissureExt_R.lvl3"
# alpha_matrix$names[alpha_matrix$names %in% c("Fimbria_L.lvl5", "Hippo_L.lvl5")] <- "Hippo_L.lvl4"
# alpha_matrix_combined <- unique(alpha_matrix)


# beta_matrix[which(my.roi$level1 == "CSF.lvl1"),]
# my.roi[which(my.roi$level2 == "CerebralCortex_R.lvl2"),]
# 
# beta_matrix[which(my.roi$level2 == "CerebralCortex_R.lvl2"),3] %>% table()
# my.roi[which(my.roi$level3 == "Frontal_R.lvl3"),2] %>% unique() %>% beta_allgroup[.,]
# 
# 
# beta_allgroup <- beta_allgroup[abs(round(beta_allgroup)) > 0, 1]
# 
# visualization.roi <- ROI.info[,-1] %>% mutate(All_brain = "all-brain")
# visualization.roi <- nz.my.roi
# q <- ncol(visualization.roi)-1
# for(j in 1:nrow(visualization.roi)){
#   string <- ""
#   for(k in (q-1):1){
#     if(!is.na(visualization.roi[j,k]) && visualization.roi[j,k] != visualization.roi[j,k+1]){
#       string <- paste(string, visualization.roi[j,k], sep = "/")
#     }
#   } 
#   # visualization.roi$pathString[j] <- gsub('.lvl[0-9]', '', string)
#   visualization.roi$pathString[j] <- string
# }
# 
# tree.structure <- as.Node(visualization.roi)
# 
# 
# gg <- as.igraph(tree.structure, directed = F, direction = "climb")
# c_scale <- colorRamp(c('red', 'white', 'blue'))
# beta_trans <- (beta[names(V(gg)),1]/0.155 + 1)/2
# V(gg)$color <- apply(c_scale(beta_trans), 1, function(x) rgb(x[1]/255,x[2]/255,x[3]/255) )
# V(gg)$label <- ifelse(abs(beta[names(V(gg)),1])>0.01, gsub('.lvl[0-9]', '', names(V(gg))), NA)
# plot(gg, vertex.size = 5)
# 
# print(tree.structure)
# 
# plot(tree.structure)
# diagonalNetwork(ToListExplicit(tree.structure, unname = TRUE), fontSize = 10, opacity = 0.9, width = 800, height=2000)
# 
# alpha_table <- data.frame(roi=ROI.info$level5[abs(alpha_disease) > 0.01], alpha=round(alpha_disease[abs(alpha_disease) > 0.01],2))
# xtable::xtable(cbind(alpha_table[1:6,], alpha_table[7:12,]))
# 
# beta_disease_coef <- cbind(thresholding(beta_disease[my.roi[,1],], 0.01), thresholding(beta_disease[my.roi[,2],], 0.01),
#                            thresholding(beta_disease[my.roi[,3],], 0.01), thresholding(beta_disease[my.roi[,4],], 0.01),
#                            thresholding(beta_disease[my.roi[,5],], 0.01), thresholding(beta_disease[my.roi[,6],], 0.01))
# roi_disease <- ROI.info[rowSums(abs(beta_disease_coef)) > 0, ]
# beta_disease_coef <- beta_disease_coef[rowSums(abs(beta_disease_coef)) > 0, ]
# c_scale <- colorRamp(c('blue', 'white', 'red'))
# beta_color <- (beta_disease_coef/0.1 + 1)/2
# data.frame(roi_disease[,4],c_scale(beta_color[,3]))
# data.frame(roi_disease[,3],c_scale(beta_color[,2]))
# data.frame(roi_disease[,2],c_scale(beta_color[,1]))
# range(beta_disease_coef)
# 
# cbind(beta1[abs(beta1$beta) > 0,1],c_scale((beta1[abs(beta1$beta) > 0,1]/0.15 + 1)/2))
# cbind(beta2[abs(beta2$beta) > 0,1],c_scale((beta2[abs(beta2$beta) > 0,1]/0.15 + 1)/2))
# cbind(beta3[abs(beta3$beta) > 0,1],c_scale((beta3[abs(beta3$beta) > 0,1]/0.15 + 1)/2))
# cbind(beta4[abs(beta4$beta) > 0,1],c_scale((beta4[abs(beta4$beta) > 0,1]/0.15 + 1)/2))
# (alpha_disease[abs(alpha_disease) > 0.01]/0.15 + 1)/2
# cbind(alpha_disease[abs(alpha_disease) > 0.01], c_scale((alpha_disease[abs(alpha_disease) > 0.01]/0.2 + 1)/2))
# 
# 
