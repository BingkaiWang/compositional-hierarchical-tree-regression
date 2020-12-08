rm(list=ls())
library(tidyverse)
library(xtable)
setwd("~/Dropbox/research/graphical-model/hierarchical-lasso/compositional-hierarchical-tree-regression/simulation/")

methods <- c("CTASSO-AIC", "CTASSO-BIC", "TASSO-AIC", "TASSO-BIC", "LASSO-AIC", "LASSO-BIC")
metrics <- c("Tuning parameter", "sensitivity", "specifity", "SSE", "df")

sim1_sl_mean <- readRDS("sim1_sl.rds") %>% map(~(apply(.,1,mean) %>% matrix(nrow = 5, ncol = 6, dimnames = list(metrics, methods)) %>% t))
sim1_st_mean <- readRDS("sim1_st.rds") %>% map(~(apply(.,1,mean) %>% matrix(nrow = 5, ncol = 6, dimnames = list(metrics, methods)) %>% t))
sim2_sl_mean <- readRDS("sim2_sl.rds") %>% map(~(apply(.[,-185],1,mean) %>% matrix(nrow = 5, ncol = 6, dimnames = list(metrics, methods)) %>% t))
sim2_st_mean <- readRDS("sim2_st.rds") %>% map(~(apply(.[,-c(31, 73)],1,mean) %>% matrix(nrow = 5, ncol = 6, dimnames = list(metrics, methods)) %>% t))
t2_mean <- rbind(sim1_sl_mean[[2]], sim1_st_mean[[2]], sim2_sl_mean[[2]], sim2_st_mean[[2]])[, c(2,3,4,1)] %>% round(2)
sim1_sl_sd <- readRDS("sim1_sl.rds") %>% map(~(apply(.,1,sd) %>% matrix(nrow = 5, ncol = 6, dimnames = list(metrics, methods)) %>% t))
sim1_st_sd <- readRDS("sim1_st.rds") %>% map(~(apply(.,1,sd) %>% matrix(nrow = 5, ncol = 6, dimnames = list(metrics, methods)) %>% t))
sim2_sl_sd <- readRDS("sim2_sl.rds") %>% map(~(apply(.[,-185],1,sd) %>% matrix(nrow = 5, ncol = 6, dimnames = list(metrics, methods)) %>% t))
sim2_st_sd <- readRDS("sim2_st.rds") %>% map(~(apply(.[,-c(31, 73)],1,sd) %>% matrix(nrow = 5, ncol = 6, dimnames = list(metrics, methods)) %>% t))
t2_sd <- rbind(sim1_sl_sd[[2]], sim1_st_sd[[2]], sim2_sl_sd[[2]], sim2_st_sd[[2]])[, c(2,3,4,1)] %>% round(2)
t2 <- data.frame(method = rownames(t2_mean),
                 tuning = rep(c("AIC", "BIC"), 12),
                 sensitivity = paste0(t2_mean[,"sensitivity"], "(", t2_sd[,"sensitivity"], ")"),
                 specifity = paste0(t2_mean[,"specifity"], "(", t2_sd[,"specifity"], ")"),
                 SSE = paste0(t2_mean[,"SSE"], "(", t2_sd[,"SSE"], ")"),
                 parameter =  paste0(t2_mean[,"Tuning parameter"], "(", t2_sd[,"Tuning parameter"], ")"))
xtable(t2)
