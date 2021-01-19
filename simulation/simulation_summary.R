rm(list=ls())
library(tidyverse)
library(xtable)
setwd("~/Dropbox/research/graphical-model/hierarchical-lasso/compositional-hierarchical-tree-regression/simulation/")

# summarizing all simulation results into tables. 
# For some settings, we remove the results of a few simulated data sets (less than 1%), since these data sets make the leaf varibles not linearly independent and the results are outliers.
methods <- c("CTASSO-AIC", "CTASSO-BIC", "TASSO-AIC", "TASSO-BIC", "LASSO-AIC", "LASSO-BIC")
metrics <- c("Tuning parameter", "sensitivity", "specifity", "SSE", "df")

# var(beta X) = var(epsilon)
sim1_sl_mean <- readRDS("sim1_sl.rds") %>% map(~(apply(.,1,mean) %>% matrix(nrow = 5, ncol = 6, dimnames = list(metrics, methods)) %>% t))
sim1_st_mean <- readRDS("sim1_st.rds") %>% map(~(apply(.,1,mean) %>% matrix(nrow = 5, ncol = 6, dimnames = list(metrics, methods)) %>% t))
sim2_sl_mean <- readRDS("sim2_sl.rds") %>% map(~(apply(.[,-c(37,548,780,967)],1,mean) %>% matrix(nrow = 5, ncol = 6, dimnames = list(metrics, methods)) %>% t))
sim2_st_mean <- readRDS("sim2_st.rds") %>% map(~(apply(.[,-c(127, 61, 247, 294, 134)],1,mean) %>% matrix(nrow = 5, ncol = 6, dimnames = list(metrics, methods)) %>% t))
t2_mean <- rbind(sim1_sl_mean[[2]], sim1_st_mean[[2]], sim2_sl_mean[[2]], sim2_st_mean[[2]])[, c(2,3,4,1)] %>% round(2)
sim1_sl_sd <- readRDS("sim1_sl.rds") %>% map(~(apply(.,1,sd) %>% matrix(nrow = 5, ncol = 6, dimnames = list(metrics, methods)) %>% t))
sim1_st_sd <- readRDS("sim1_st.rds") %>% map(~(apply(.,1,sd) %>% matrix(nrow = 5, ncol = 6, dimnames = list(metrics, methods)) %>% t))
sim2_sl_sd <- readRDS("sim2_sl.rds") %>% map(~(apply(.[,-c(37,548,780,967)],1,sd) %>% matrix(nrow = 5, ncol = 6, dimnames = list(metrics, methods)) %>% t))
sim2_st_sd <- readRDS("sim2_st.rds") %>% map(~(apply(.[,-c(127, 61, 247, 294, 134)],1,sd) %>% matrix(nrow = 5, ncol = 6, dimnames = list(metrics, methods)) %>% t))
t2_sd <- rbind(sim1_sl_sd[[2]], sim1_st_sd[[2]], sim2_sl_sd[[2]], sim2_st_sd[[2]])[, c(2,3,4,1)] %>% round(2)
t2 <- data.frame(method = rownames(t2_mean),
                 tuning = rep(c("AIC", "BIC"), 12),
                 sensitivity = paste0(t2_mean[,"sensitivity"], "(", t2_sd[,"sensitivity"], ")"),
                 specifity = paste0(t2_mean[,"specifity"], "(", t2_sd[,"specifity"], ")"),
                 SSE = paste0(t2_mean[,"SSE"], "(", t2_sd[,"SSE"], ")"),
                 parameter =  paste0(t2_mean[,"Tuning parameter"], "(", t2_sd[,"Tuning parameter"], ")"))
t2
xtable(t2)

# 10var(beta X) = var(epsilon)
sim1_sl_mean <- readRDS("sim1_sl.rds") %>% map(~(apply(.,1,mean) %>% matrix(nrow = 5, ncol = 6, dimnames = list(metrics, methods)) %>% t))
sim1_st_mean <- readRDS("sim1_st.rds") %>% map(~(apply(.[,-c(436, 115, 412, 625, 882, 187, 419)],1,mean) %>% matrix(nrow = 5, ncol = 6, dimnames = list(metrics, methods)) %>% t))
sim2_sl_mean <- readRDS("sim2_sl.rds") %>% map(~(apply(.[,-c(875,58, 609,967, 626, 70, 930)],1,mean) %>% matrix(nrow = 5, ncol = 6, dimnames = list(metrics, methods)) %>% t))
sim2_st_mean <- readRDS("sim2_st.rds") %>% map(~(apply(.[,-c(5,21,39)],1,mean) %>% matrix(nrow = 5, ncol = 6, dimnames = list(metrics, methods)) %>% t))
t3_mean <- rbind(sim1_sl_mean[[3]], sim1_st_mean[[3]], sim2_sl_mean[[3]], sim2_st_mean[[3]])[, c(2,3,4,1)] %>% round(2)
sim1_sl_sd <- readRDS("sim1_sl.rds") %>% map(~(apply(.,1,sd) %>% matrix(nrow = 5, ncol = 6, dimnames = list(metrics, methods)) %>% t))
sim1_st_sd <- readRDS("sim1_st.rds") %>% map(~(apply(.[,-c(436, 115, 412, 625, 882, 187, 419)],1,sd) %>% matrix(nrow = 5, ncol = 6, dimnames = list(metrics, methods)) %>% t))
sim2_sl_sd <- readRDS("sim2_sl.rds") %>% map(~(apply(.[,-c(875,58, 609,967, 626, 70, 930)],1,sd) %>% matrix(nrow = 5, ncol = 6, dimnames = list(metrics, methods)) %>% t))
sim2_st_sd <- readRDS("sim2_st.rds") %>% map(~(apply(.[,-c(5,21,39)],1,sd) %>% matrix(nrow = 5, ncol = 6, dimnames = list(metrics, methods)) %>% t))
t3_sd <- rbind(sim1_sl_sd[[3]], sim1_st_sd[[3]], sim2_sl_sd[[3]], sim2_st_sd[[3]])[, c(2,3,4,1)] %>% round(2)
t3 <- data.frame(method = rownames(t3_mean),
                 tuning = rep(c("AIC", "BIC"), 12),
                 sensitivity = paste0(t3_mean[,"sensitivity"], "(", t3_sd[,"sensitivity"], ")"),
                 specifity = paste0(t3_mean[,"specifity"], "(", t3_sd[,"specifity"], ")"),
                 SSE = paste0(t3_mean[,"SSE"], "(", t3_sd[,"SSE"], ")"),
                 parameter =  paste0(t3_mean[,"Tuning parameter"], "(", t3_sd[,"Tuning parameter"], ")"))
t3
xtable(t3)
