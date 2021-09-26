rm(list=ls())
setwd("~/Dropbox\ (Personal)/research/graphical-model/hierarchical-lasso/compositional-hierarchical-tree-regression/data-analysis/")
set.seed(123)
library(tidyverse)
library(genlasso)
library(foreach)
library(doParallel)
library(tictoc)
cl <- makeCluster(4)
registerDoParallel(cl)
packages <- c("tidyverse", "genlasso", "MASS")
load("ADNI/MRICloud/200925/Data_vol_v1_Type1.RData")
MRIcloud_table <- readxl::read_xlsx("ADNI/MRICloud_table.xlsx")
my.roi <- MRIcloud_table %>% 
  as.data.frame %>%
  filter(Exclude == 0) %>% 
  .[,2:6] %>%
  mutate(All_brain = "all-brain")

# construct_D is a function to calculate D(eta) for CT-LASSO of the main paper. 
# construct_D_pairwise is a function to calculate D(eta) for CT-LASSO-p of the main paper. 
# construct_TASSO is a function to construct the regularization function of TASSO.
# calculate_beta is a function to calculate beta based on alpha.
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
construct_D_pairwise <- function(my.roi, eta, scale = F, sd = NULL){
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
  D <- matrix(0, nrow = q*(q-1)/2, ncol = q)
  colnames(D) <- my.roi[,1]
  iter <- 1
  for(k in 1:(ncol(my.roi)-1)){
    for(i in unique(my.roi[,k+1])){
      index <- my.roi[my.roi[,k+1] == i,k]
      if(length(unique(index)) > 1){
        for(j in 1:(length(unique(index))-1)){
          leafindex_j <- my.roi[my.roi[,k] == unique(index)[j],1]
          for(jj in (j+1):length(unique(index))){
            leafindex_jj <- my.roi[my.roi[,k] == unique(index)[jj],1]
            D[iter, leafindex_j] <- t(weight.mat[my.roi[,k] == unique(index)[j],k])
            D[iter, leafindex_jj] <- -t(weight.mat[my.roi[,k] == unique(index)[jj],k])
            iter <- iter + 1
          }
        }
      }
    }
  }
  D <- D[1:(iter-1),]
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
thresholding <- function(x, threshold) {ifelse(abs(x) > threshold, x, 0)}


# data preprocessing =============
# my.roi <- ROI.info[,-c(1,7)] %>% mutate(All_brain = "all-brain")
my.roi[73:74,3] <- c("Amyg_L", "Amyg_R") # correct ROI names such that no region has more than 1 parent
# my.roi[89:90, 2] <- c("BasalForebrainBG_L", "BasalForebrainBG_R")
my.roi[107:108, 2] <- c("PVA_inferior_L", "PVA_inferior_R")
for(j in 1:5){my.roi[,j] <- paste0(my.roi[,j], ".lvl", 6-j)} # Since different regions may have the same name, this is to tell them apart.
for(j in 5:2){
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
X_complete <- dat[[6]][nonmissing_indi, MRIcloud_table$Exclude == 0]
X_complete <- X_complete / (rowSums(X_complete) %*% t(rep(1, ncol(X_complete))))
n <- length(Y_complete)
saveRDS(list(X_complete, my.roi), "../sim2-data.rds")
saveRDS(my.roi, "../tree-structure.rds")

# All group analysis CTASSO ==============
tictoc::tic()
param_grid <- expand.grid(eta = seq(0, 1, by = 0.001))
tuning_param_all <- foreach(j = 1:nrow(param_grid), .combine = cbind, .packages = packages) %dopar% {
  D_mat <- construct_D(my.roi, param_grid[j,1])
  genlasso.fit <- genlasso(y = Y_complete, X = X_complete, D = D_mat)
  df <- genlasso.fit$df
  temp <- rep(NA, 600)
  temp[1:length(df)] <- n * log(colSums((Y_complete %*% t(rep(1,length(df))) - genlasso.fit$fit)^2)) + log(n) * df
  temp
}
tuning_param <-tuning_param_all %>% apply(.,2,min, na.rm = T) %>% which.min %>% param_grid[.,] %>% unlist()
D_mat <- construct_D(my.roi, tuning_param[1])
genlasso.fit <- genlasso(y = Y_complete, X = X_complete, D = D_mat)
df <- genlasso.fit$df
IC <- n * log((colSums((Y_complete %*% t(rep(1,length(df))) - genlasso.fit$fit)^2))) + log(n) * df
alpha_allgroup <- coef(genlasso.fit, lambda = genlasso.fit$lambda[which.min(IC)])$beta
alpha_allgroup <- alpha_allgroup - mean(alpha_allgroup)
beta_allgroup <- calculate_beta(my.roi, alpha_allgroup)
tictoc::toc()
# saveRDS(list(tuning_param_all, alpha_allgroup, beta_allgroup), file = "ADNI-CTASSO-BIC-allgroup.rds")
# write.table(cbind(thresholding(alpha_allgroup, 0.01), MRIcloud_table$index[MRIcloud_table$Exclude==0]), "alpha-allgroup.txt", col.names = F)
# write.table(cbind(thresholding(beta_allgroup, 0.01), 1:length(beta_allgroup)), "beta-allgroup.txt", col.names = F)

# All group analysis TASSO ============
D_mat <- construct_TASSO(my.roi)
genlasso.fit <- genlasso(y = Y_complete, X = X_complete, D = D_mat)
df <- genlasso.fit$df
IC <- n * log((colSums((scale(Y_complete, scale = F) %*% t(rep(1,length(df))) - genlasso.fit$fit)^2))) + log(n) * df
alpha_TASSO <- coef(genlasso.fit, lambda = genlasso.fit$lambda[which.min(IC)])$beta
beta_TASSO <- calculate_beta(my.roi, alpha_TASSO)
sum(round(alpha_TASSO - mean(alpha_TASSO)) > 0)

# all group analysis CLASSO ==========
genlasso.fit <- genlasso(y = Y_complete, X = X_complete, D = construct_D(my.roi, 1))
df <- genlasso.fit$df
IC <- n * log((colSums((scale(Y_complete, scale = F) %*% t(rep(1,length(df))) - genlasso.fit$fit)^2))) + log(n) * df
alpha_LASSO <- coef(genlasso.fit, lambda = genlasso.fit$lambda[which.min(IC)])$beta
beta_LASSO <- calculate_beta(my.roi, alpha_LASSO)
data.frame(names = my.roi[,1], alpha = alpha_LASSO - mean(alpha_LASSO))
sum(round(alpha_LASSO - mean(alpha_LASSO)) > 0)

# all group analysis TFL =========
genlasso.fit <- genlasso(y = Y_complete, X = X_complete, D = construct_D(my.roi, 0))
df <- genlasso.fit$df
IC <- n * log((colSums((scale(Y_complete, scale = F) %*% t(rep(1,length(df))) - genlasso.fit$fit)^2))) + log(n) * df
alpha_TFL <- coef(genlasso.fit, lambda = genlasso.fit$lambda[which.min(IC)])$beta
beta_TFL <- calculate_beta(my.roi, alpha_TFL)
data.frame(names = my.roi[,1], alpha = alpha_TFL - mean(alpha_TFL))
sum(round(alpha_TFL - mean(alpha_TFL)) > 0)

# all group analysis CT-LASSO-pairwise ==========
tictoc::tic()
param_grid <- expand.grid(eta = seq(0, 1, by = 0.001))
tuning_param_all <- foreach(j = 1:nrow(param_grid), .combine = cbind, .packages = packages) %dopar% {
  D_mat <- construct_D_pairwise(my.roi, param_grid[j,1])
  genlasso.fit <- genlasso(y = Y_complete, X = X_complete, D = D_mat)
  df <- genlasso.fit$df
  temp <- rep(NA, 600)
  temp[1:length(df)] <- n * log(colSums((Y_complete %*% t(rep(1,length(df))) - genlasso.fit$fit)^2)) + log(n) * df
  temp
}
tuning_param <-tuning_param_all %>% apply(.,2,min, na.rm = T) %>% which.min %>% param_grid[.,] %>% unlist()
D_mat <- construct_D_pairwise(my.roi, tuning_param[1])
genlasso.fit <- genlasso(y = Y_complete, X = X_complete, D = D_mat)
df <- genlasso.fit$df
IC <- n * log((colSums((Y_complete %*% t(rep(1,length(df))) - genlasso.fit$fit)^2))) + log(n) * df
alpha_pairwise <- coef(genlasso.fit, lambda = genlasso.fit$lambda[which.min(IC)])$beta
alpha_pairwise <- alpha_pairwise - mean(alpha_pairwise)
beta_pairwise <- calculate_beta(my.roi, alpha_pairwise)
tictoc::toc()

# Comparison of CTASSO, TASSO, CLASSO, TFL, CLASSO-pairwise ==========
alpha_comparison <- data.frame(names = my.roi[,1], 
                               alpha_CTASSO = round(alpha_allgroup - mean(alpha_allgroup),2),
                               alpha_TASSO = round(alpha_TASSO - mean(alpha_TASSO),2),
                               alpha_LASSO = round(alpha_LASSO - mean(alpha_LASSO),2),
                               alpha_TFL = round(alpha_TFL - mean(alpha_TFL),2)
                               # alpha_pairwise = round(alpha_pairwise - mean(alpha_pairwise),2)
                               )
beta_comparison <- data.frame(beta_CTASSO = round(beta_allgroup - mean(beta_allgroup),2),
                             beta_TASSO = round(beta_TASSO - mean(beta_TASSO),2),
                             beta_LASSO = round(beta_LASSO - mean(beta_LASSO),2),
                             beta_TFL = round(beta_TFL - mean(beta_TFL),2)
                             )
colSums(abs(beta_comparison)>0)
(abs(alpha_comparison[,2:5]) > 0) %>% apply(2,sum)
(abs(beta_comparison) > 0) %>% apply(2,sum)
(abs(alpha_comparison[,2:5]) > 0) %>% apply(1,sum) %>% table
(abs(beta_comparison) > 0) %>% apply(1,sum) %>% table


# data.frame(names = my.roi[,1], alpha = alpha_allgroup * apply(X_complete, 2, sd) * 1000)


# plot(lm(Y_complete~., data.frame(cbind(Y_complete, X_complete[,"SylParieSul_L"]))))
# 
# ggplot(data.frame(cbind(memory_score = Y_complete, lingual_gyrus = X_complete[,"LG_L"]))) + 
#   geom_point(aes(x = lingual_gyrus, y = memory_score)) +
#   geom_smooth(method = lm, formula = y~x,aes(x = lingual_gyrus, y = memory_score)) +
#   annotate("text", x = 0.0105, y = 3, label = "Slope = -122.00 (P-value 0.002)", size = 5)
# 
# ggplot(data.frame(cbind(memory_score = Y_complete, cuneus = X_complete[,"Cu_L"]))) + 
#   geom_point(aes(x = cuneus, y = memory_score)) +
#   geom_smooth(method = lm, formula = y~x,aes(x = cuneus, y = memory_score) ) +
#   annotate("text", x = 0.007, y = 3, label = "Slope = -63.88 (P-value 0.161)", size = 5)
# 
# 
# ggplot(data.frame(cbind(memory_score = Y_complete, SylParieSul = X_complete[,"SylParieSul_L"]))) + 
#   geom_point(aes(x = SylParieSul, y = memory_score)) +
#   geom_smooth(method = lm, formula = y~x,aes(x = SylParieSul, y = memory_score) ) + 
#   annotate("text", x = 0.00075, y = 3, label = "Slope = -362.91 (P-value 0.074)", size = 5)

