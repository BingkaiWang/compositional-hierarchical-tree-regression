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

# CTASSO with AIC and scaling X ------------
tuning_param_all <- foreach(j = 1:nrow(param_grid), .combine = cbind, .packages = packages) %dopar% {
  D_mat <- construct_D(my.roi, param_grid[j,1], scale = T, sd = apply(X,2,sd) * 1e5)
  genlasso.fit <- genlasso(y = Y.centered, X = scale(X), D = D_mat, eps = param_grid[j,2])
  df <- genlasso.fit$df - n
  n * log(colSums((Y.centered %*% t(rep(1,length(df))) - genlasso.fit$fit)^2)) + 2 * df
}
tuning_param <-tuning_param_all %>% apply(.,2,min) %>% which.min %>% param_grid[.,] %>% unlist()
D_mat <- construct_D(my.roi, tuning_param[1], scale = T, sd = apply(X,2,sd) * 1e5)
genlasso.fit <- genlasso(y = Y.centered, X = scale(X), D = D_mat, eps = tuning_param[2])
df <- genlasso.fit$df - n
IC <- n * log((colSums((Y.centered %*% t(rep(1,length(df))) - genlasso.fit$fit)^2))) + 2 * df
alpha <- coef(genlasso.fit, lambda = genlasso.fit$lambda[which.min(IC)])$beta
beta <- calculate_beta(my.roi, alpha)

print(" CTASSO with AIC and scaling X ")
print(tuning_param)
saveRDS(list(tuning_param_all, alpha, beta), file = "CTASSO-AIC.rds")

# CTASSO with BIC and scaling X ------------
tuning_param_all <- foreach(j = 1:nrow(param_grid), .combine = cbind, .packages = packages) %dopar% {
  D_mat <- construct_D(my.roi, param_grid[j,1], scale = T, sd = apply(X,2,sd) * 1e5)
  genlasso.fit <- genlasso(y = Y.centered, X = scale(X), D = D_mat, eps = param_grid[j,2])
  df <- genlasso.fit$df - n
  n * log(colSums((Y.centered %*% t(rep(1,length(df))) - genlasso.fit$fit)^2)) + log(n) * df
}
tuning_param <-tuning_param_all %>% apply(.,2,min) %>% which.min %>% param_grid[.,] %>% unlist()
D_mat <- construct_D(my.roi, tuning_param[1], scale = T, sd = apply(X,2,sd) * 1e5)
genlasso.fit <- genlasso(y = Y.centered, X = scale(X), D = D_mat, eps = tuning_param[2])
df <- genlasso.fit$df - n
IC <- n * log((colSums((Y.centered %*% t(rep(1,length(df))) - genlasso.fit$fit)^2))) + log(n) * df
alpha <- coef(genlasso.fit, lambda = genlasso.fit$lambda[which.min(IC)])$beta
beta <- calculate_beta(my.roi, alpha)

print(" CTASSO with BIC and scaling X ")
print(tuning_param)
saveRDS(list(tuning_param_all, alpha, beta), file = "CTASSO-BIC.rds")


# TASSO with AIC and scaling X ------------
tuning_param_all <- map(1:nrow(param_grid_tasso), function(j) {
  genlasso.fit <- genlasso(y = Y.centered, X = scale(X), D = construct_TASSO(my.roi), eps = param_grid_tasso[j,1])
  df <- genlasso.fit$df - n
  n * log(colSums((Y.centered %*% t(rep(1,length(df))) - genlasso.fit$fit)^2)) + 2 * df
})
tuning_param <-tuning_param_all %>% map_dbl(min) %>% which.min %>% param_grid_tasso[.,] %>% unlist()
genlasso.fit <- genlasso(y = Y.centered, X = scale(X), D = construct_TASSO(my.roi), eps = tuning_param)
df <- genlasso.fit$df - n
IC <- n * log((colSums((Y.centered %*% t(rep(1,length(df))) - genlasso.fit$fit)^2))) + 2 * df
alpha <- coef(genlasso.fit, lambda = genlasso.fit$lambda[which.min(IC)])$beta
beta <- calculate_beta(my.roi, alpha)

print(" TASSO with AIC and scaling X ")
saveRDS(list(tuning_param_all, alpha, beta), file = "TASSO-AIC.rds")

# TASSO with BIC and scaling X ------------
tuning_param_all <- map(1:nrow(param_grid_tasso), function(j) {
  genlasso.fit <- genlasso(y = Y.centered, X = scale(X), D = construct_TASSO(my.roi), eps = param_grid_tasso[j,1])
  df <- genlasso.fit$df - n
  n * log(colSums((Y.centered %*% t(rep(1,length(df))) - genlasso.fit$fit)^2)) + log(n) * df
})
tuning_param <-tuning_param_all %>% map_dbl(min) %>% which.min %>% param_grid_tasso[.,] %>% unlist()
genlasso.fit <- genlasso(y = Y.centered, X = scale(X), D = construct_TASSO(my.roi), eps = tuning_param)
df <- genlasso.fit$df - n
IC <- n * log((colSums((Y.centered %*% t(rep(1,length(df))) - genlasso.fit$fit)^2))) + log(n) * df
alpha <- coef(genlasso.fit, lambda = genlasso.fit$lambda[which.min(IC)])$beta
beta <- calculate_beta(my.roi, alpha)

print(" TASSO with BIC and scaling X ")
print(tuning_param)
saveRDS(list(tuning_param_all, alpha, beta), file = "TASSO-BIC.rds")



# LASSO with AIC and scaling X ------------
tuning_param_all <- map(1:nrow(param_grid_lasso), function(j) {
  genlasso.fit <- genlasso(y = Y.centered, X = scale(X), D = diag(nrow(my.roi)), eps = param_grid_lasso[j,1])
  df <- genlasso.fit$df - n
  n * log(colSums((Y.centered %*% t(rep(1,length(df))) - genlasso.fit$fit)^2)) + 2 * df
})
tuning_param <-tuning_param_all %>% map_dbl(min) %>% which.min %>% param_grid_lasso[.,] %>% unlist()
genlasso.fit <- genlasso(y = Y.centered, X = scale(X), D = diag(nrow(my.roi)), eps = tuning_param)
df <- genlasso.fit$df - n
IC <- n * log((colSums((Y.centered %*% t(rep(1,length(df))) - genlasso.fit$fit)^2))) + 2 * df
alpha <- coef(genlasso.fit, lambda = genlasso.fit$lambda[which.min(IC)])$beta
beta <- calculate_beta(my.roi, alpha)

print(" LASSO with AIC and scaling X ")
print(tuning_param)
saveRDS(list(tuning_param_all, alpha, beta), file = "LASSO-AIC.rds")

# LASSO with BIC and scaling X ------------
tuning_param_all <- map(1:nrow(param_grid_lasso), function(j) {
  genlasso.fit <- genlasso(y = Y.centered, X = scale(X), D = diag(nrow(my.roi)), eps = param_grid_lasso[j,1])
  df <- genlasso.fit$df - n
  n * log(colSums((Y.centered %*% t(rep(1,length(df))) - genlasso.fit$fit)^2)) + log(n) * df
})
tuning_param <-tuning_param_all %>% map_dbl(min) %>% which.min %>% param_grid_lasso[.,] %>% unlist()
genlasso.fit <- genlasso(y = Y.centered, X = scale(X), D = diag(nrow(my.roi)), eps = tuning_param)
df <- genlasso.fit$df - n
IC <- n * log((colSums((Y.centered %*% t(rep(1,length(df))) - genlasso.fit$fit)^2))) + log(n) * df
alpha <- coef(genlasso.fit, lambda = genlasso.fit$lambda[which.min(IC)])$beta
beta <- calculate_beta(my.roi, alpha)

print(" LASSO with BIC and scaling X ")
print(tuning_param)
saveRDS(list(tuning_param_all, alpha, beta), file = "LASSO-BIC.rds")



# CTASSO with 10-fold CV and scaling X --------
n_fold <- 10
test_size <- floor(n/n_fold)
param_grid <- expand.grid(lambda = exp(seq(-1,7, by = 0.05)), eta = seq(0, 1, by = 0.05), gamma = 1e-4)
random_order <- sample(1:n, replace = F)
random_partition <- map(1:n_fold, ~random_order[(.-1)*test_size + 1:test_size])
random_partition[[n_fold]] <- random_order[((n_fold-1)*test_size + 1):n]
X.scaled <- scale(X)
tuning_param_all <- foreach(j = 1:nrow(param_grid), .combine = cbind, .packages = packages) %dopar% {
  map_dbl(1:n_fold, function(k){
    Y_train <- scale(Y[-random_partition[[k]]], scale = F)
    X_train <- X.scaled[-random_partition[[k]],]
    Y_test <- scale(Y[random_partition[[k]]], scale = F)
    X_test <- X.scaled[random_partition[[k]], ]
    D_mat <- construct_D(my.roi, param_grid[j,2], scale = T, sd = apply(X,2,sd) * 1e5)
    genlasso.fit <- genlasso(y = Y_train, X = X_train, D = D_mat, eps = param_grid[j,3])
    if(param_grid[j,1] > min(genlasso.fit$lambda)){
      genlasso.pred <- predict(genlasso.fit, lambda = param_grid[j,1], Xnew = X_test)
    }else{
      write.table(NA, file=paste0(round(param_grid[j,1],3),"-replaced-by-", round(min(genlasso.fit$lambda),3),".txt"))
      genlasso.pred <- predict(genlasso.fit, lambda = max(genlasso.fit$lambda), Xnew = X_test)
    }
    sum((Y_test - genlasso.pred$fit)^2)
  }) %>% mean
}
tuning_param <-tuning_param_all %>% which.min %>% param_grid[.,] %>% unlist()

print("CTASSO with 10-fold CV and scaling X")
print(tuning_param)

D_mat <- construct_D(my.roi, tuning_param[2], scale = T, sd = apply(X,2,sd) * 1e5)
genlasso.fit <- genlasso(y = Y.centered, X = X.scaled, D = D_mat, eps = tuning_param[3])
alpha <- coef(genlasso.fit, lambda = tuning_param[1])$beta
beta <- calculate_beta(my.roi, alpha)

saveRDS(list(tuning_param_all, alpha, beta), file = "CTASSO-CV.rds")

# visualization -----
# coef.roi4 <- mutate(my.roi, coef = readRDS("data_analysis4.rds")[[2]]) %>% mutate(size = colMeans(data.level5))
# coef.roi2 <- mutate(my.roi, coef = readRDS("data_analysis3.rds")[[2]]) %>% mutate(size = colMeans(data.level5))
# 
# dd <- data.frame(leaf_nodes = my.roi$`level5-type1`, 
#                  AIC = as.vector(round(readRDS("data_analysis3.rds")[[2]], 3)),
#                  CV = as.vector(round(readRDS("data_analysis4.rds")[[2]], 3)))
# dd[dd$CV != 0, c(1,3) ] %>% xtable()
# 
# treemap(coef.roi4,
#         index = colnames(my.roi)[5:1],
#         vSize = "size", vColor = "coef", type = "value")
# dd <- coef.roi4[abs(coef.roi4$coef) > 0.01,c(1,7)]
# dd$coef <- round(dd$coef, 2)

# 
# # coef.roi1 <- mutate(my.roi, coef = round(coef.fit1$beta,4)) %>% mutate(size = colMeans(data.level5))
# treemap(coef.roi1,
#         index = colnames(my.roi)[5:1],
#         vSize = "size", vColor = "coef", type = "value")

# # Model fit with AIC without scaling X ------------
# param_grid <- expand.grid(eta = seq(0, 1, by = 0.05), gamma = c(1e-4, 1e-2))
# tuning_param_all <- foreach(j = 1:nrow(param_grid), .combine = cbind, .packages = packages) %dopar% {
#   genlasso.fit <- genlasso(y = Y.centered, X = scale(X, scale = F), D = construct_D(my.roi, param_grid[j,1]), eps = param_grid[j,2])
#   df <- genlasso.fit$df - n
#   n * log(colSums((Y.centered %*% t(rep(1,length(df))) - genlasso.fit$fit)^2)) + 2 * df
# }
# 
# tuning_param <-tuning_param_all %>% apply(.,2,min) %>% which.min %>% param_grid[.,] %>% unlist()
# genlasso.fit1 <- genlasso(y = Y.centered, X = scale(X, scale = F), D = construct_D(my.roi, tuning_param[1]), eps = tuning_param[2])
# df <- genlasso.fit1$df - n
# IC <- n * log((colSums((Y.centered %*% t(rep(1,length(df))) - genlasso.fit1$fit)^2))) + 2 * df
# alpha1 <- coef(genlasso.fit1, lambda = genlasso.fit1$lambda[which.min(IC)])$beta
# beta1 <- calculate_beta(my.roi, alpha1)
# 
# print("Model fit with AIC without scaling X")
# print(tuning_param)
# saveRDS(list(tuning_param_all, alpha1, beta1), file = "data_analysis1.rds")
# 
# 
# # Model fit with 10-fold CV without scaling X --------
# n_fold <- 10
# test_size <- floor(n/n_fold)
# param_grid <- expand.grid(eta = seq(0, 1, by = 0.05), gamma = 1e-4, lambda = exp(seq(-6,1, by = 0.02)))
# random_order <- sample(1:n, replace = F)
# random_partition <- map(1:n_fold, ~random_order[(.-1)*test_size + 1:test_size])
# random_partition[[n_fold]] <- random_order[((n_fold-1)*test_size + 1):n]
# tuning_param_all <- foreach(j = 1:nrow(param_grid), .combine = cbind, .packages = packages) %dopar% {
#   map_dbl(1:n_fold, function(k){
#     Y_train <- scale(Y[-random_partition[[k]]], scale = F)
#     X_train <- scale(X[-random_partition[[k]],], scale = F)
#     Y_test <- scale(Y[random_partition[[k]]], scale = F)
#     X_test <- scale(X[random_partition[[k]], ], scale = F)
#     genlasso.fit <- genlasso(y = Y_train, X = X_train, D = construct_D(my.roi, param_grid[j,1]), eps = param_grid[j,2])
#     if(param_grid[j,3] > min(genlasso.fit$lambda)){
#       genlasso.pred <- predict(genlasso.fit, lambda = param_grid[j,3], Xnew = X_test)
#     }else{
#       write.table(NA, file=paste0(round(param_grid[j,3],3),"-replaced-by-", round(min(genlasso.fit$lambda),3),".txt"))
#       genlasso.pred <- predict(genlasso.fit, lambda = min(genlasso.fit$lambda), Xnew = X_test)
#     }
#     sum((Y_test - genlasso.pred$fit)^2)
#   }) %>% mean
# }
# tuning_param <-tuning_param_all %>% which.min %>% param_grid[.,] %>% unlist()
# 
# print("Model fit with 10-fold CV with scaling X")
# print(tuning_param)
# 
# genlasso.fit2 <- genlasso(y = Y.centered, X = scale(X, scale = F), D = construct_D(my.roi, tuning_param[1]), eps = tuning_param[2])
# alpha2 <- coef(genlasso.fit2, lambda = tuning_param[3])$beta
# beta2 <- calculate_beta(my.roi, alpha2)
# 
# saveRDS(list(tuning_param_all, alpha2, beta2), file = "data_analysis2.rds")
# 
# # # 
# # # 
# # # 
# # # fit2 <- genlasso(y = Y.centered, X = scale(X, scale = F), D = construct_D(my.roi,0.5), eps = 1)
# # # coef.fit2 <- coef(fit2, lambda=sqrt(length(Y)*log(ncol(X))))
# # # plot(coef.fit2$beta)
# # # sum(abs(coef.fit2$beta) > 0.0001)
# # # plot(fit2)
# # # coef.roi2 <- mutate(my.roi, coef = round(coef.fit2$beta,4)) %>% mutate(size = colMeans(data.level5))
# # # treemap(coef.roi2, 
# # #         index = colnames(my.roi)[5:1],
# # #         vSize = "size", vColor = "coef", type = "value")
# # # 
# # # 
# # # fit2 <- genlasso(y = Y.centered, X = scale(X, scale = F), D = construct_D(my.roi,0.999), eps = 1)
# # # coef.fit2 <- coef(fit2, lambda=sqrt(length(Y)*log(ncol(X))))
# # # plot(coef.fit2$beta)
# # # sum(abs(coef.fit2$beta) > 0.0001)
# # # plot(fit2)
# # # coef.roi2 <- mutate(my.roi, coef = round(coef.fit2$beta,4)) %>% mutate(size = colMeans(data.level5))
# # # treemap(coef.roi2, 
# # #         index = colnames(my.roi)[5:1],
# # #         vSize = "size", vColor = "coef", type = "value")
# # # 
# # # 
# # # fit2 <- genlasso(y = Y.centered, X = scale(X), D = diag(nrow(my.roi)), eps = 1)
# # # 
# # # 
# # # fit3 <- genlasso(y = scale(Y), X = scale(X), D = construct_D(my.roi,0.5), eps = 1)
# # # coef.fit3 <- coef(fit3, lambda=sqrt(length(Y)*log(ncol(X))))
# # # plot(coef.fit3$beta)
# # # sum(abs(coef.fit3$beta) > 0.0001)
# # # plot(fit3)
# # # coef.roi2 <- mutate(my.roi, coef = round(coef.fit3$beta,4)) %>% mutate(size = colMeans(data.level5))
# # # treemap(coef.roi2,
# # #         index = colnames(my.roi)[5:1],
# # #         vSize = "size", vColor = "coef", type = "value")
# # # 
# # # 
# # validating compositional structure ----------------
# map_dbl(unique(my.roi[,5]), function(j){
#   d1 <- data.level1[,j]
#   d2 <- data.level2[, colnames(data.level2) %in% my.roi[my.roi[,5] == j,4]]
#   abs(rowSums(d2) - d1) %>% mean()
# }) %>% table
# map_dbl(unique(my.roi[,4]), function(j){
#   d2 <- data.level2[,j]
#   d3 <- data.level3[, colnames(data.level3) %in% my.roi[my.roi[,4] == j,3]]
#   if(is.null(dim(d3))){
#     abs(d3 - d2) %>% mean()
#   }else{
#     abs(rowSums(d3) - d2) %>% mean()
#   }
# }) %>% round(4) %>% table
# map_dbl(unique(my.roi[,3]), function(j){
#   d3 <- data.level3[,j]
#   d4 <- data.level4[, colnames(data.level4) %in% my.roi[my.roi[,3] == j,2]]
#   if(is.null(dim(d4))){
#     abs(d4 - d3) %>% mean()
#   }else{
#     abs(rowSums(d4) - d3) %>% mean()
#   }
# }) %>% round(4) %>% table
# map_dbl(unique(my.roi[,2]), function(j){
#   d4 <- data.level4[,j]
#   d5 <- data.level5[, colnames(data.level5) %in% my.roi[my.roi[,2] == j,1]]
#   if(is.null(dim(d5))){
#     mean(abs(d5 - d4))
#   }else{
#     mean(abs(rowSums(d5) - d4))
#   }
# }) %>% round(4) %>% table()
# 
# 
# data.level1 <- dat.level.type1[["level1"]][,unique(my.roi[,5])]
# data.level2 <- dat.level.type1[["level2"]][,unique(my.roi[,4])]
# data.level5 <- dat.level.type1[["level5"]][,my.roi[,1]]
# missing.indi <- rowSums(leaf.data) %>% is.na | is.na(Y)
# leaf.data <- leaf.data[!missing.indi,]
# for(i in 1:nrow(leaf.data)){leaf.data[i,] <- leaf.data[i,]/sum(leaf.data[i,])}
# Y <- Y[!missing.indi]
# 
# # method 2: generalized LASSO without composition
# 
# # # method 3: LASSO with composition
# # map_dbl(1:ncol(leaf.data), function(i){sum(leaf.data[,i]<1e-9)})
# # fit.3 <- cv.glmnet(x = log(leaf.data), y = Y)
# 
# 
# 
# 
# 
# # visualization
# library(ggraph)
# my.roi <- ROI.info[,2:6] %>% filter(! `level1-type1` %in% c(NA, "Myelencephalon")) %>% distinct()
# write.csv(my.roi, "my.roi.csv")
# roi.edges <- rbind(cbind("1", my.roi[,5]) %>% data.frame() %>% distinct(),
#                    cbind( my.roi[,5], my.roi[,4]) %>% data.frame() %>% distinct(),
#                    cbind( my.roi[,4], my.roi[,3]) %>% data.frame() %>% distinct(),
#                    cbind( my.roi[,3], my.roi[,2]) %>% data.frame() %>% distinct(),
#                    cbind( my.roi[,2], my.roi[,1]) %>% data.frame() %>% distinct())
# roi.edges <- roi.edges[which(map_dbl(1:nrow(roi.edges), ~roi.edges[.,1] == roi.edges[.,2]) == 0), ] %>% data.frame() %>% distinct() 
# # which(map_dbl(unique(roi.edges[,2]), ~nrow(roi.edges[roi.edges[,2] == .,])) == 2)
# ggraph(graph_from_data_frame( roi.edges ), layout = 'dendrogram', circular = FALSE) + 
#   geom_edge_diagonal() +
#   geom_node_point() +
#   theme_void()
# 
# 
# 
# ########
# my.roi[my.roi[,3] == "BasalGang_L", ]
# 
# (data.level4[,"Caud_L"] + data.level4[,"Caudate_tail_L"] + data.level4[,"GP_L"] + data.level4[,"Put_L"] + data.level4[,"BasalForebrain_L"] - data.level3[,"BasalGang_L"]) %>% plot()
# 
# (data.level3[,"BasalGang_R"] + data.level4[,"Amyg_R"] - data.level2[,"CerebralNucli_R"]) %>% plot()
# 
# data.level3[,"BasalGang_L"] %>% plot
# 
# (data.level4[,"Caud_L"] + data.level4[,"Caudate_tail_L"] + data.level4[,"GP_L"] + data.level4[,"Put_L"] + data.level5[,"NucAccumbens_L"] - data.level3[,"BasalGang_L"]) %>% plot()
# 
# (data.level4[,"PVA_posterior_L"] - data.level5[,"PVWp_L"] - data.level5[,"PVWl_L"]) %>% plot
# 
# 
# mm <- colnames(data.level4) %in% my.roi[my.roi[,3] == "Limbic_L",2]
# rowSums(data.level4[,mm]) - data.level3[,"Limbic_L"]
# 
# plot(data.level5[,"Hippo_R"] + data.level5[,"Fimbria_R"] - data.level4[,"Hippo_R"])
# 
# plot(data.level5[,"Hippo_R"] - data.level4[,"Hippo_R"])
# 