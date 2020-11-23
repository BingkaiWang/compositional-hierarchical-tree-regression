rm(list=ls())
setwd("~/Dropbox/research/graphical-model/hierarchical-lasso/compositional-hierarchical-tree-regression/data-analysis/")
set.seed(123)
library(tidyverse)
library(genlasso)
library(foreach)
library(doParallel)
cl <- makeCluster(2)
registerDoParallel(cl)
packages <- c("tidyverse", "genlasso", "MASS")
load("ADNI/MRICloud/200925/Data_vol_v1_Type1.RData")
construct_D <- function(my.roi, eta, scale = F, sd = NULL){
  q <- nrow(my.roi)
  weight.mat <- matrix(1, nrow = q, ncol = ncol(my.roi)-1)
  for(k in 2:ncol(weight.mat)){
    for(i in unique(my.roi[,k])){
      weight.mat[my.roi[,k] == i, k] <- 1/length(unique(my.roi[my.roi[,k] == i, k-1]))
    }
  }
  for(i in 1:nrow(weight.mat)){weight.mat[i,] <- cumprod(weight.mat[i,])}
  if(scale == T){
    weight.mat <- diag(sd) %*% weight.mat
  }
  D <- matrix(0, nrow = q-1, ncol = q)
  colnames(D) <- my.roi[,1]
  iter <- 1
  for(k in 1:(ncol(my.roi)-1)){
    for(i in unique(my.roi[,k+1])){
      index <- my.roi[my.roi[,k+1] == i,k]
      if(length(unique(index)) > 1){
        for(j in 1:(length(unique(index))-1)){
          leafindex_j <- my.roi[my.roi[,k] == unique(index)[j],1]
          leafindex_jp1 <- my.roi[my.roi[,k] == unique(index)[j+1],1]
          D[iter, leafindex_j] <- t(weight.mat[my.roi[,k] == unique(index)[j],k])
          D[iter, leafindex_jp1] <- -t(weight.mat[my.roi[,k] == unique(index)[j+1],k])
          iter <- iter + 1
        }
      }
    }
  }
  if(eta == 0){
    D_mat <- D
  } else if(eta == 1){
    D_mat <- (diag(q) - matrix(1/q, nrow = q, ncol = q))
  } else {
    D_mat <- rbind((diag(q) - matrix(1/q, nrow = q, ncol = q)) * eta, D * (1-eta))
  }
  return(D_mat)
}
calculate_beta <- function(my.roi, alpha){
  my.roi <- my.roi[,-ncol(my.roi)]
  p <- length(unique(unlist(my.roi)))
  p.leaf <- length(alpha)
  C <- matrix(0, nrow = p, ncol = p)
  colnames(C) <- unique(unlist(my.roi))
  u <- C[1,]
  for(i in 1:p.leaf){
    C[i, unique(as.character(my.roi[i,]))] <- 1
  }
  mark <- p.leaf + 1
  for(i in 2:ncol(my.roi)){
    index <- unique(my.roi[,i])
    for(j in index){
      constraint.var <- as.character(my.roi[my.roi[,i] == j,i-1])
      if(length(unique(constraint.var)) > 1){
        C[mark, unique(constraint.var)] <- 1
        mark <- mark + 1
      }
    }
  }
  u[as.character(my.roi[,ncol(my.roi)])] <- 1
  tilde_alpha <- c(alpha, rep(0, p-p.leaf))
  tilde_1 <- c(rep(1, p.leaf), rep(0, p-p.leaf))
  inv_C <- solve(C)
  c <- as.numeric(u %*% inv_C %*% tilde_alpha/(u %*% inv_C %*% tilde_1))
  beta <- inv_C %*% (tilde_alpha - c * tilde_1)
  return(beta)
}
thresholding <- function(x, threshold) {ifelse(abs(x) > threshold, x, 0)}

my.roi <- ROI.info[,-1] %>% mutate(All_brain = "all-brain")
my.roi[73:74,3] <- c("LimbicCN_L", "LimbicCN_R") # correct ROI names such that no region has more than 1 parent
my.roi[89:90, 2] <- c("BasalForebrainBG_L", "BasalForebrainBG_R")
my.roi[129:130, 2] <- c("PVA_posteriorIWM_L", "PVA_posteriorIWM_R")
for(j in 1:6){my.roi[,j] <- paste0(my.roi[,j], ".lvl", 6-j)} # Since different regions may have the same name, this is to tell them apart.
for(j in 6:2){
  index <- unique(my.roi[,j])
  for(k in index){
    if(length(unique(my.roi[my.roi[,j] == k,j-1])) == 1){my.roi[my.roi[,j] == k,j-1] <- k}
  }
}
Y <- dat.cog[,"ADNI_MEM"]
Diagnosis <- dat.demo$Dx
nonmissing_indi <- complete.cases(cbind(Y, Diagnosis, dat[[6]]))
Y_complete <- Y[nonmissing_indi]
Diagnosis_complete <- Diagnosis[nonmissing_indi]
X_complete <- dat[[6]][nonmissing_indi, ]
X_complete <- X_complete / (rowSums(X_complete) %*% t(rep(1, ncol(X_complete))))
param_grid <- expand.grid(eta = seq(0, 1, by = 0.05), gamma = c(1e-4, 1e-2))
n <- length(Y_complete)
# saveRDS(list(X_complete, my.roi), "../simulation/sim2-data.rds")

# All group analysis
tuning_param_all <- foreach(j = 1:nrow(param_grid), .combine = cbind, .packages = packages) %dopar% {
  D_mat <- construct_D(my.roi, param_grid[j,1], scale = T, sd = apply(X_complete,2,sd) * 1e5)
  genlasso.fit <- genlasso(y = scale(Y_complete), X = scale(X_complete), D = D_mat, eps = param_grid[j,2])
  df <- genlasso.fit$df
  n * log(colSums((scale(Y_complete) %*% t(rep(1,length(df))) - genlasso.fit$fit)^2)) + log(n) * df
}
tuning_param <-tuning_param_all %>% apply(.,2,min) %>% which.min %>% param_grid[.,] %>% unlist()
D_mat <- construct_D(my.roi, tuning_param[1], scale = T, sd = apply(X_complete,2,sd) * 1e5)
genlasso.fit <- genlasso(y = scale(Y_complete), X = scale(X_complete), D = D_mat, eps = tuning_param[2])
df <- genlasso.fit$df
IC <- n * log((colSums((scale(Y_complete) %*% t(rep(1,length(df))) - genlasso.fit$fit)^2))) + log(n) * df
alpha_allgroup <- coef(genlasso.fit, lambda = genlasso.fit$lambda[which.min(IC)])$beta
beta_allgroup <- calculate_beta(my.roi, alpha_allgroup)
saveRDS(list(tuning_param_all, alpha_allgroup, beta_allgroup), file = "ADNI-CTASSO-BIC-allgroup.rds")
write.table(cbind(thresholding(alpha_allgroup, 0.01), ROI.info$Index), "alpha-allgroup.txt", col.names = F)
# write.table(cbind(thresholding(beta_allgroup, 0.01), 1:length(beta_allgroup)), "beta-allgroup.txt", col.names = F)

# AD + MCI group analysis
Y_disease <- Y_complete[Diagnosis_complete != "CN"]
X_disease <- X_complete[Diagnosis_complete != "CN",]
X_disease <- X_disease / (rowSums(X_disease) %*% t(rep(1, ncol(X_disease))))
n <- length(Y_disease)
tuning_param_all <- foreach(j = 1:nrow(param_grid), .combine = cbind, .packages = packages) %dopar% {
  D_mat <- construct_D(my.roi, param_grid[j,1], scale = T, sd = apply(X_disease,2,sd) * 1e5)
  genlasso.fit <- genlasso(y = scale(Y_disease), X = scale(X_disease), D = D_mat, eps = param_grid[j,2])
  df <- genlasso.fit$df
  n * log(colSums((scale(Y_disease) %*% t(rep(1,length(df))) - genlasso.fit$fit)^2)) + log(n) * df
}
tuning_param <-tuning_param_all %>% apply(.,2,min) %>% which.min %>% param_grid[.,] %>% unlist()
D_mat <- construct_D(my.roi, tuning_param[1], scale = T, sd = apply(X_disease,2,sd) * 1e5)
genlasso.fit <- genlasso(y = scale(Y_disease), X = scale(X_disease), D = D_mat, eps = tuning_param[2])
df <- genlasso.fit$df
IC <- n * log((colSums((scale(Y_disease) %*% t(rep(1,length(df))) - genlasso.fit$fit)^2))) + log(n) * df
alpha_disease <- coef(genlasso.fit, lambda = genlasso.fit$lambda[which.min(IC)])$beta
beta_disease <- calculate_beta(my.roi, alpha_disease)
saveRDS(list(tuning_param_all, alpha_disease, beta_disease), file = "ADNI-CTASSO-BIC-disease.rds")
write.table(cbind(thresholding(alpha_disease, 0.01), ROI.info$Index), "alpha-disease.txt", col.names = F)
# write.table(cbind(thresholding(beta_disease, 0.01), 1:length(beta_disease)), "beta-disease.txt", col.names = F)
beta_disease_coef <- cbind(thresholding(beta_disease[my.roi[,1],], 0.01), thresholding(beta_disease[my.roi[,2],], 0.01),
                           thresholding(beta_disease[my.roi[,3],], 0.01), thresholding(beta_disease[my.roi[,4],], 0.01),
                           thresholding(beta_disease[my.roi[,5],], 0.01), thresholding(beta_disease[my.roi[,6],], 0.01))
beta1 <- data.frame(beta=ifelse(my.roi$level3 == "LateralVentricle_L.lvl3", thresholding(rowSums(beta_disease_coef), 0.01), 0), ROI.info$Index)
write.table(beta1, "beta1.txt", col.names = F)
beta2 <- data.frame(beta=ifelse(my.roi$level4 == "BasalForebrain_L.lvl2", thresholding(rowSums(beta_disease_coef), 0.01), 0), ROI.info$Index)
write.table(beta2, "beta2.txt", col.names = F)
beta3 <- data.frame(beta=ifelse(my.roi$level3 %in% c("Temporal_L.lvl3", "Occipital_L.lvl3"), thresholding(rowSums(beta_disease_coef), 0.01), 0), ROI.info$Index)
write.table(beta3, "beta3.txt", col.names = F)
beta4 <- data.frame(beta=ifelse(my.roi$level3 %in% c("Temporal_R.lvl3", "Occipital_R.lvl3"), thresholding(rowSums(beta_disease_coef), 0.01), 0), ROI.info$Index)
write.table(beta4, "beta4.txt", col.names = F)

s# data.frame(roi = my.roi[,1], coef = alpha) %>% filter(abs(alpha) > 0.001)
# data.frame(roi = unique(unlist(my.roi[,-7])), coef = beta) %>% filter(abs(beta) > 0.01)
# 
# x5 <- thresholding(beta[my.roi[,1],], 0.01)
# x4 <- round(beta[unique(my.roi[,2]),], 2)
# 
# beta_coef <- cbind(thresholding(beta[my.roi[,1],], 0.01), thresholding(beta[my.roi[,2],], 0.01), 
#                    thresholding(beta[my.roi[,3],], 0.01), thresholding(beta[my.roi[,4],], 0.01), 
#                    thresholding(beta[my.roi[,5],], 0.01), thresholding(beta[my.roi[,6],], 0.01))
# 
# nz.my.roi <- my.roi[apply(abs(beta_coef), 1, sum) > 0, ]
# nz.beta_coef <- beta_coef[apply(abs(beta_coef), 1, sum) > 0, ]
# for(i in 1:(ncol(nz.my.roi)-1)){
#   for(j in 1:nrow(nz.my.roi)){
#     if(sum(abs(nz.beta_coef[j,1:i])) == 0){nz.my.roi[j,1:i] <- NA}
#   }
# }

