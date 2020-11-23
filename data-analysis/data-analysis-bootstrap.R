rm(list=ls())
set.seed(123)
# setwd("~/Dropbox/research/graphical-model/hierarchical-lasso/compositional-hierarchical-tree-regression/data-analysis/")
setwd("~/graphical_model/hierarchical-lasso")
load("MRICloud/200408/Data.RData")
library(tidyverse)
library(glmnet)
library(genlasso)
library(treemap)
library(MASS)
library(foreach)
library(doParallel)
cl <- makeCluster(16)
registerDoParallel(cl)
packages <- c("tidyverse", "genlasso", "MASS")
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
  D <- vector("list", ncol(my.roi)-1)
  for(k in 1:length(D)){
    D[[k]] <- bdiag(0)
    for(i in unique(my.roi[,k+1])){
      index <- my.roi[my.roi[,k+1] == i,k]
      if(length(unique(index)) > 1){
        d <- matrix(0, nrow = length(unique(index))-1, ncol = length(index))
        for(j in 1:nrow(d)){
          d[j,index == unique(index)[j]] <- t(weight.mat[my.roi[,k] == unique(index)[j],k])
          d[j,index == unique(index)[j+1]] <- -t(weight.mat[my.roi[,k] == unique(index)[j+1],k])
        }
        D[[k]] <- bdiag(D[[k]], d)
      } else {
        D[[k]] <- bdiag(D[[k]], matrix(0, nrow = 1, ncol = length(index)))
      }
    }
    D[[k]] <- D[[k]][rowSums(abs(D[[k]]))>0, -1]
  }
  D_mat <- (diag(q) - matrix(1/q, nrow = q, ncol = q)) * eta
  for(k in 1:length(D)){
    D_mat <- rbind(D_mat, D[[k]] * (1-eta))
  }
  D_mat <- D_mat[rowSums(abs(D_mat))>0, ]
  return(D_mat)
}
construct_TASSO <- function(my.roi){
  q <- nrow(my.roi)
  p <- length(unique(unlist(my.roi)))
  P_star <- diag(q) - matrix(1/q, nrow = q, ncol = q)
  P_2 <- matrix(0, nrow = p-q, ncol = q)
  mark <- 1
  for(k in 2:ncol(my.roi)){
    for(i in unique(my.roi[,k])){
      index <- my.roi[my.roi[,k] == i, 1:(k-1), drop = F] %>% apply(., 2, function(x){length(unique(x))})
      if(min(index) > 1){
        P_2[mark, which(my.roi[,k] == i)] <- 1
        mark <- mark + 1
      }
    }
  }
  if(mark != p-q+1){
    stop("An error has occurred!")
  } else{
    return(rbind(0.5 * P_star, 0.5 * P_2 %*% P_star))
  }
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

# Type 1 analysis ==================
# Preprocessing ------------
# define outcome variable
Y <- dat.demo[,"PicVocab_Unadj"]
# remove Myelencephalon from brain regions & remove duplicates & reorganize data to achieve composition
# The following brain regions have multiple parents: PVA_posterior_L, PVA_posterior_R, BasalForebrain_L, BasalForebrain_R, Limbic_L and Limbic_R
# Add the following regions at Level 5 to achieve compoistion: BFRemainder_L, BFRemainder_R, HippoRemainder_L, HippoRemainder_R, Limbic_L and Limbic_R
my.roi <- ROI.info[,2:6] %>% filter(! `level1-type1` %in% c(NA, "Myelencephalon")) %>% distinct() %>%
  arrange(`level1-type1`, `level2-type1`, `level3-type1`, `level4-type1`, `level5-type1`)
my.roi[148,2] <- "PVWl_L"
my.roi[170,2] <- "PVWp_L"
my.roi[248,2] <- "PVWl_R"
my.roi[270,2] <- "PVWp_R"
my.roi[109,2] <- "NucAccumbens_L"
my.roi[209,2] <- "NucAccumbens_R"
my.roi[114,3] <- "Amyg_L"
my.roi[214,3] <- "Amyg_R"
dat.level.type1[["level4"]] <- cbind(dat.level.type1[["level4"]], dat.level.type1[["level5"]][, c("PVWl_L", "PVWp_L", "PVWl_R", "PVWp_R", "NucAccumbens_L", "NucAccumbens_R")])
dat.level.type1[["level4"]][,"BasalForebrain_L"] <- dat.level.type1[["level3"]][,"BasalForebrain_L"]
dat.level.type1[["level4"]][,"BasalForebrain_R"] <- dat.level.type1[["level3"]][,"BasalForebrain_R"]
dat.level.type1[["level3"]] <- cbind(dat.level.type1[["level3"]], dat.level.type1[["level4"]][,c("Amyg_L", "Amyg_R")])
dat.level.type1[["level3"]][,"Limbic_L"] <- dat.level.type1[["level3"]][,"Limbic_L"] - dat.level.type1[["level4"]][,"Amyg_L"]
dat.level.type1[["level3"]][,"Limbic_R"] <- dat.level.type1[["level3"]][,"Limbic_R"] - dat.level.type1[["level4"]][,"Amyg_R"]
my.roi <- rbind(my.roi,
                c("BFRemainder_R","BasalForebrain_R", "BasalForebrain_R", "BasalForebrain_R","Diencephalon_R"),
                c("BFRemainder_L","BasalForebrain_L", "BasalForebrain_L", "BasalForebrain_L","Diencephalon_L"),
                c("HippoRemainder_L", "Hippo_L", "Limbic_L", "CerebralCortex_L", "Telencephalon_L"),
                c("HippoRemainder_R", "Hippo_R", "Limbic_R", "CerebralCortex_R", "Telencephalon_R"),
                c("Limbic_L", "Limbic_L", "Limbic_L", "CerebralCortex_L", "Telencephalon_L"),
                c("Limbic_R", "Limbic_R", "Limbic_R", "CerebralCortex_R", "Telencephalon_R"))
dat.level.type1[["level5"]] <- cbind(dat.level.type1[["level5"]],
                                     BFRemainder_L = dat.level.type1[["level4"]][,"BasalForebrain_L"] - rowSums(dat.level.type1[["level5"]][,c("BasalForebrain_L", "Cl_L", "HypoThalamus_L", "Mammillary_L")]),
                                     BFRemainder_R = dat.level.type1[["level4"]][,"BasalForebrain_R"] - rowSums(dat.level.type1[["level5"]][,c("BasalForebrain_R", "Cl_R", "HypoThalamus_R", "Mammillary_R")]),
                                     HippoRemainder_L = dat.level.type1[["level4"]][,"Hippo_L"] - rowSums(dat.level.type1[["level5"]][,c("Fimbria_L", "Hippo_L")]),
                                     HippoRemainder_R = dat.level.type1[["level4"]][,"Hippo_R"] - rowSums(dat.level.type1[["level5"]][,c("Fimbria_R", "Hippo_R")]),
                                     Limbic_L = dat.level.type1[["level4"]][,"Limbic_L"] - rowSums(dat.level.type1[["level5"]][,c("ENT_L", "PHG_L")]),
                                     Limbic_R = dat.level.type1[["level4"]][,"Limbic_R"] - rowSums(dat.level.type1[["level5"]][,c("ENT_R", "PHG_R")])
)
data.level5 <- dat.level.type1[["level5"]][,my.roi[,1]]
complete.indi <- complete.cases(data.level5) & !is.na(Y)
data.level1 <- dat.level.type1[["level1"]][complete.indi,unique(my.roi[,5])]
data.level2 <- dat.level.type1[["level2"]][complete.indi,unique(my.roi[,4])]
data.level3 <- dat.level.type1[["level3"]][complete.indi,unique(my.roi[,3])]
data.level4 <- dat.level.type1[["level4"]][complete.indi,unique(my.roi[,2])]
data.level5 <- data.level5[complete.indi,]
my.roi <- cbind(my.roi, `level0-type1`=rep("0", nrow(my.roi)))
for(j in 1:5){my.roi[,j] <- paste0(my.roi[,j], ".lvl", 6-j)} # Since different regions may have the same name, this is to tell them apart.
for(j in 5:2){
  index <- unique(my.roi[,j])
  for(k in index){
    if(length(unique(my.roi[my.roi[,j] == k,j-1])) == 1){my.roi[my.roi[,j] == k,j-1] <- k}
  }
}
Y <- Y[complete.indi]
X <- data.level5 / (rowSums(data.level5) %*% t(rep(1, ncol(data.level5))))
Y.centered <- scale(Y, scale = F)
n <- length(Y)
param_grid <- expand.grid(eta = seq(0, 1, by = 0.05), gamma = c(1e-4, 1e-2))
param_grid_lasso <- expand.grid(gamma = c(1e-4, 1e-2))
param_grid_tasso <- expand.grid(gamma = c(1e-4, 1e-2))


# CTASSO with AIC and scaling X ------------
n_bootstrap <- 1000
CTASSO_AIC_boot <- foreach(j = 1:n_bootstrap, .combine = cbind, .packages = packages, .errorhandling = "remove") %dopar% {
  index <- sample(1:n, size = 0.8 * n, replace = F)
  Y.bootstrap <- Y.centered[index]
  X.bootstrap <- X[index,]
  tuning_param_all <- map(1:nrow(param_grid), function(j) {
    D_mat <- construct_D(my.roi, param_grid[j,1], scale = T, sd = apply(X.bootstrap,2,sd) * 1e5)
    genlasso.fit <- genlasso(y = Y.bootstrap, X = scale(X.bootstrap), D = D_mat, eps = param_grid[j,2])
    df <- genlasso.fit$df - n
    n * log(colSums((Y.centered %*% t(rep(1,length(df))) - genlasso.fit$fit)^2)) + 2 * df
  })
  tuning_param <-tuning_param_all %>% map_dbl(min) %>% which.min %>% param_grid[.,] %>% unlist()
  D_mat <- construct_D(my.roi, tuning_param[1], scale = T, sd = apply(X.bootstrap,2,sd) * 1e5)
  genlasso.fit <- genlasso(y = Y.bootstrap, X = scale(X.bootstrap), D = D_mat, eps = tuning_param[2])
  df <- genlasso.fit$df - n
  IC <- n * log((colSums((Y.bootstrap %*% t(rep(1,length(df))) - genlasso.fit$fit)^2))) + 2 * df
  alpha <- coef(genlasso.fit, lambda = genlasso.fit$lambda[which.min(IC)])$beta
  alpha
}
saveRDS(CTASSO_AIC_boot, file = "CTASSO-AIC-bootstrap.rds")

# CTASSO with BIC and scaling X ------------
n_bootstrap <- 1000
CTASSO_BIC_boot <- foreach(j = 1:n_bootstrap, .combine = cbind, .packages = packages, .errorhandling = "remove") %dopar% {
  index <- sample(1:n, size = 0.8 * n, replace = F)
  Y.bootstrap <- Y.centered[index]
  X.bootstrap <- X[index,]
  tuning_param_all <- map(1:nrow(param_grid), function(j) {
    D_mat <- construct_D(my.roi, param_grid[j,1], scale = T, sd = apply(X.bootstrap,2,sd) * 1e5)
    genlasso.fit <- genlasso(y = Y.bootstrap, X = scale(X.bootstrap), D = D_mat, eps = param_grid[j,2])
    df <- genlasso.fit$df - n
    n * log(colSums((Y.centered %*% t(rep(1,length(df))) - genlasso.fit$fit)^2)) + log(n) * df
  })
  tuning_param <-tuning_param_all %>% map_dbl(min) %>% which.min %>% param_grid[.,] %>% unlist()
  D_mat <- construct_D(my.roi, tuning_param[1], scale = T, sd = apply(X.bootstrap,2,sd) * 1e5)
  genlasso.fit <- genlasso(y = Y.bootstrap, X = scale(X.bootstrap), D = D_mat, eps = tuning_param[2])
  df <- genlasso.fit$df - n
  IC <- n * log((colSums((Y.bootstrap %*% t(rep(1,length(df))) - genlasso.fit$fit)^2))) + log(n) * df
  alpha <- coef(genlasso.fit, lambda = genlasso.fit$lambda[which.min(IC)])$beta
  alpha
}
saveRDS(CTASSO_BIC_boot, file = "CTASSO-BIC-bootstrap.rds")


