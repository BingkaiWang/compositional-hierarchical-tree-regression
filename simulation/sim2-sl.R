rm(list=ls())
set.seed(123)
# setwd("~/Dropbox/research/graphical-model/hierarchical-lasso/compositional-hierarchical-tree-regression/simulation/")
setwd("~/graphical_model/hierarchical-lasso")
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
    D_mat <- rbind(0.5 * P_star, 0.5 * P_2 %*% P_star)
    return(D_mat[-nrow(D_mat),])
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


# data setup =====================
n_sim <- 1000
X_complete <- readRDS("sim2-data.rds")[[1]]
sim.tree <- readRDS("sim2-data.rds")[[2]]
n <- nrow(X_complete)
p.leaf <- nrow(sim.tree)
param_grid <- expand.grid(eta = seq(0.0, 1, by = 0.05))
param_grid_lasso <- expand.grid(gamma = c(1e-4, 1e-2))
param_grid_tasso <- expand.grid(gamma = c(1e-4, 1e-2))
noise_sl <- sd(3 * X_complete[,1] - 2 * X_complete[,3] - X_complete[,5]) * sqrt(c(0.1, 1, 10)) # sd(beta X):sd(eplsilon) = 5, 1, 0.2

# sparse tree: Y = 3 * SFG_L - 2 * SFG_PFC_L - SFG_pole_L
beta_true <- c(3, 0, -2, 0, -1, 0, rep(0,363))
sim2_sl <- vector("list", length(noise_sl))
for(t in 1:length(sim2_sl)){
  sim2_sl[[t]] <- foreach(j = 1:n_sim, .combine = cbind, .packages = packages) %dopar% {
    X.leaf <- X_complete[sample(1:n, n, replace = T),]
    X.leaf.centered <- scale(X.leaf, scale = F)
    Y <- 3 + 3 * X.leaf.centered[,1] - 2 * X.leaf.centered[,3] - X.leaf.centered[,5] + rnorm(n, sd = noise_sl[t])
    Y.centered <- scale(Y, scale = F)
    result <- rep(NA, length = 30)
    
    # CTASSO with AIC
    result[1:5] <- tryCatch(
      {tunning_param <- map_dbl(1:nrow(param_grid), function(j){
        sim.fit <- genlasso(y = Y, X = X.leaf, D = construct_D(sim.tree, param_grid[j,1]))
        df <- sim.fit$df
        IC <- n * log(colSums((Y %*% t(rep(1,length(df))) - sim.fit$fit)^2)) + 2 * df
        min(IC)
      }) %>% which.min %>% param_grid[.,] %>% unlist()
      sim.fit <- genlasso(y = Y, X = X.leaf, D = construct_D(sim.tree, tunning_param[1]))
      df <- sim.fit$df
      IC <- n * log(colSums((Y %*% t(rep(1,length(df))) - sim.fit$fit)^2)) + 2 * df
      alpha <- coef(sim.fit, lambda = sim.fit$lambda[which.min(IC)])$beta
      beta <- calculate_beta(sim.tree, alpha)
      c(tunning_param[1], sensitivity = mean(abs(beta[which(abs(beta_true)>0)]) > 0.01), 
        specifity = mean(abs(beta[which(abs(beta_true)==0)]) < 0.01), MSE = sum((beta - beta_true)^2), ENP = sum(abs(beta) > 0.01))},
      error = function(err) {print(err); rep(NA,5)}
    )
    
    # CTASSO with BIC
    result[6:10] <- tryCatch(
      {tunning_param <- map_dbl(1:nrow(param_grid), function(j){
        sim.fit <- genlasso(y = Y, X = X.leaf, D = construct_D(sim.tree, param_grid[j,1]))
        df <- sim.fit$df
        IC <- n * log(colSums((Y %*% t(rep(1,length(df))) - sim.fit$fit)^2)) + log(n) * df
        min(IC)
      }) %>% which.min %>% param_grid[.,] %>% unlist()
      sim.fit <- genlasso(y = Y, X = X.leaf, D = construct_D(sim.tree, tunning_param[1]))
      df <- sim.fit$df
      IC <- n * log(colSums((Y %*% t(rep(1,length(df))) - sim.fit$fit)^2)) + log(n) * df
      alpha <- coef(sim.fit, lambda = sim.fit$lambda[which.min(IC)])$beta
      beta <- calculate_beta(sim.tree, alpha)
      c(tunning_param[1], sensitivity = mean(abs(beta[which(abs(beta_true)>0)]) > 0.01),
        specifity = mean(abs(beta[which(abs(beta_true)==0)]) < 0.01), MSE = sum((beta - beta_true)^2), ENP = sum(abs(beta) > 0.01))},
      error = function(err) {print(err); rep(NA,5)}
    )
    
    # TASSO with AIC
    result[11:15] <- tryCatch(
      {
        tunning_param <- map_dbl(1:nrow(param_grid_tasso), function(j){
          sim.fit <- genlasso(y = Y, X = X.leaf.centered, D = construct_TASSO(sim.tree), eps = param_grid_tasso[j,1])
          df <- sim.fit$df
          IC <- n * log(colSums((Y.centered %*% t(rep(1,length(df))) - sim.fit$fit)^2)) + 2 * df
          min(IC)
        }) %>% which.min %>% param_grid_tasso[.,] %>% unlist()
        sim.fit <- genlasso(y = Y, X = X.leaf.centered, D = construct_TASSO(sim.tree), eps = tunning_param)
        df <- sim.fit$df
        IC <- n * log(colSums((Y.centered %*% t(rep(1,length(df))) - sim.fit$fit)^2)) + 2 * df
        alpha <- coef(sim.fit, lambda = sim.fit$lambda[which.min(IC)])$beta
        beta <- calculate_beta(sim.tree, alpha)
        c(tunning_param, sensitivity = mean(abs(beta[which(abs(beta_true)>0)]) > 0.01), 
          specifity = mean(abs(beta[which(abs(beta_true)==0)]) < 0.01), MSE = sum((beta - beta_true)^2), ENP = sum(abs(beta) > 0.01))
      },
      error = function(err) {print(err); rep(NA,5)}
    )
    
    # TASSO with BIC
    result[16:20] <- tryCatch(
      {
        tunning_param <- map_dbl(1:nrow(param_grid_tasso), function(j){
          sim.fit <- genlasso(y = Y, X = X.leaf.centered, D = construct_TASSO(sim.tree), eps = param_grid_tasso[j,1])
          df <- sim.fit$df
          IC <- n * log(colSums((Y.centered %*% t(rep(1,length(df))) - sim.fit$fit)^2)) + log(n) * df
          min(IC)
        }) %>% which.min %>% param_grid_tasso[.,] %>% unlist()
        sim.fit <- genlasso(y = Y, X = X.leaf.centered, D = construct_TASSO(sim.tree), eps = tunning_param)
        df <- sim.fit$df
        IC <- n * log(colSums((Y.centered %*% t(rep(1,length(df))) - sim.fit$fit)^2)) + log(n) * df
        alpha <- coef(sim.fit, lambda = sim.fit$lambda[which.min(IC)])$beta
        beta <- calculate_beta(sim.tree, alpha)
        c(tunning_param, sensitivity = mean(abs(beta[which(abs(beta_true)>0)]) > 0.01), 
          specifity = mean(abs(beta[which(abs(beta_true)==0)]) < 0.01), MSE = sum((beta - beta_true)^2), ENP = sum(abs(beta) > 0.01))
      },
      error = function(err) {print(err); rep(NA,5)}
    )
    
    # LASSO with AIC
    result[21:25] <- tryCatch(
      {
        tunning_param <- map_dbl(1:nrow(param_grid_lasso), function(j){
          sim.fit <- genlasso(y = Y, X = X.leaf.centered, D = diag(p.leaf), eps = param_grid_lasso[j,1])
          df <- sim.fit$df
          IC <- n * log(colSums((Y.centered %*% t(rep(1,length(df))) - sim.fit$fit)^2)) + 2 * df
          min(IC)
        }) %>% which.min %>% param_grid_lasso[.,] %>% unlist()
        
        sim.fit <- genlasso(y = Y, X = X.leaf.centered, D = diag(p.leaf), eps = tunning_param)
        df <- sim.fit$df
        IC <- n * log(colSums((Y.centered %*% t(rep(1,length(df))) - sim.fit$fit)^2)) + 2 * df
        alpha <- coef(sim.fit, lambda = sim.fit$lambda[which.min(IC)])$beta
        beta <- calculate_beta(sim.tree, alpha)
        c(tunning_param, sensitivity = mean(abs(beta[which(abs(beta_true)>0)]) > 0.01), 
          specifity = mean(abs(beta[which(abs(beta_true)==0)]) < 0.01), MSE = sum((beta - beta_true)^2), ENP = sum(abs(beta) > 0.01))
      },
      error = function(err) {print(err); rep(NA,5)}
    )
    
    # LASSO with BIC
    result[26:30] <- tryCatch(
      {
        tunning_param <- map_dbl(1:nrow(param_grid_lasso), function(j){
          sim.fit <- genlasso(y = Y, X = X.leaf.centered, D = diag(p.leaf), eps = param_grid_lasso[j,1])
          df <- sim.fit$df
          IC <- n * log(colSums((Y.centered %*% t(rep(1,length(df))) - sim.fit$fit)^2)) + log(n) * df
          min(IC)
        }) %>% which.min %>% param_grid_lasso[.,] %>% unlist()
        
        sim.fit <- genlasso(y = Y, X = X.leaf.centered, D = diag(p.leaf), eps = tunning_param)
        df <- sim.fit$df
        IC <- n * log(colSums((Y.centered %*% t(rep(1,length(df))) - sim.fit$fit)^2)) + log(n) * df
        alpha <- coef(sim.fit, lambda = sim.fit$lambda[which.min(IC)])$beta
        beta <- calculate_beta(sim.tree, alpha)
        c(tunning_param, sensitivity = mean(abs(beta[which(abs(beta_true)>0)]) > 0.01), 
          specifity = mean(abs(beta[which(abs(beta_true)==0)]) < 0.01), MSE = sum((beta - beta_true)^2), ENP = sum(abs(beta) > 0.01))
      },
      error = function(err) {print(err); rep(NA,5)}
    )
    result
  }
}
saveRDS(object = sim2_sl, file = "sim2_sl.rds")

