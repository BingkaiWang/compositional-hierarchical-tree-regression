rm(list=ls())
library(tidyverse)
library(xtable)
setwd("~/Dropbox/research/graphical-model/hierarchical-lasso/compositional-hierarchical-tree-regression/simulation/")

methods <- c("CTASSO-AIC", "CTASSO-BIC", "TASSO-AIC", "TASSO-BIC", "LASSO-AIC", "LASSO-BIC")
metrics <- c("Tuning parameter", "sensitivity", "specifity", "MSE", "df")

sim1_sl <- readRDS("sim1_sl.rds") %>% map(~apply(.,1,mean) %>% matrix(nrow = 5, ncol = 6, dimnames = list(metrics, methods)) %>% t)
sim1_st <- readRDS("sim1_st.rds") %>% map(~apply(.,1,mean) %>% matrix(nrow = 5, ncol = 6, dimnames = list(metrics, methods)) %>% t)
sim2_sl <- readRDS("sim2_sl.rds") %>% map(~apply(.,1,mean) %>% matrix(nrow = 5, ncol = 6, dimnames = list(metrics, methods)) %>% t)
sim2_st <- readRDS("sim2_st.rds") %>% map(~apply(.,1,mean) %>% matrix(nrow = 5, ncol = 6, dimnames = list(metrics, methods)) %>% t)


xtable(rbind(sim1_sl[[2]][,-1], sim1_st[[2]][,-1], sim2_sl[[2]][,-1], sim2_st[[2]][,-1]))

