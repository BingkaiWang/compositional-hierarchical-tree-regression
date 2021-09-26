rm(list=ls())
library(tidyverse)
library(xtable)
library(cowplot)
setwd("~/Dropbox (Personal)/research/graphical-model/hierarchical-lasso/compositional-hierarchical-tree-regression/")

# summarizing all simulation results into tables. 
# For some settings, we remove the results of a few simulated data sets (less than 1%), 
# since these data sets make the leaf variables not linearly independent (our assumption violated)
# and the results are outliers.
methods <- c("CTASSO-AIC", "CTASSO-BIC", "TASSO-AIC", "TASSO-BIC", "LASSO-AIC", "LASSO-BIC", "TFL-2-AIC", "TFL-2-BIC", "CTASSOp-AIC", "CTASSOp-BIC")
metrics <- c("Tuning parameter", "sensitivity", "specifity", "l2loss", "MSE", "df")

# =====================
# var(beta X) = var(epsilon)
c1 <- c(260, 345, 413, 612, 682, 709, 750, 911)
c2 <- c(206, 231, 360, 553, 606, 620)
sim1_sl <- readRDS("sim1_sl.rds")
sim1_st <- readRDS("sim1_st.rds")
sim2_sl <- readRDS("sim2_sl.rds")
sim2_st <- readRDS("sim2_st.rds")
for(i in 1:length(sim1_sl)){
  sim1_sl[[i]][seq(4, 60, by = 6), ] <- sqrt(sim1_sl[[i]][seq(4, 60, by = 6), ])
  sim1_st[[i]][seq(4, 60, by = 6), ] <- sqrt(sim1_st[[i]][seq(4, 60, by = 6), ])
  sim2_sl[[i]][seq(4, 60, by = 6), ] <- sqrt(sim2_sl[[i]][seq(4, 60, by = 6), ])
  sim2_st[[i]][seq(4, 60, by = 6), ] <- sqrt(sim2_st[[i]][seq(4, 60, by = 6), ])
}
sim1_sl_mean <- sim1_sl %>% map(~(apply(.,1,mean) %>% matrix(nrow = 6, ncol = 10, dimnames = list(metrics, methods)) %>% t))  %>% map(~.[c(1,2,9,10,7,8,3,4,5,6),])
sim1_st_mean <- sim1_st %>% map(~(apply(.,1,mean) %>% matrix(nrow = 6, ncol = 10, dimnames = list(metrics, methods)) %>% t)) %>% map(~.[c(1,2,9,10,7,8,3,4,5,6),])
sim2_sl_mean <- sim2_sl %>% map(~(apply(.[,-c1],1,mean, na.rm = T) %>% matrix(nrow = 6, ncol = 10, dimnames = list(metrics, methods)) %>% t)) %>% map(~.[c(1,2,9,10,7,8,3,4,5,6),])
sim2_st_mean <- sim2_st %>% map(~(apply(.[,-c2],1,mean, na.rm = T) %>% matrix(nrow = 6, ncol = 10, dimnames = list(metrics, methods)) %>% t)) %>% map(~.[c(1,2,9,10,7,8,3,4,5,6),])
t2_mean <- rbind(sim1_sl_mean[[1]], sim1_st_mean[[1]], sim2_sl_mean[[1]], sim2_st_mean[[1]])[, c(2,3,4,1)] %>% round(2)
sim1_sl_sd <- sim1_sl %>% map(~(apply(.,1,sd) %>% matrix(nrow = 6, ncol = 10, dimnames = list(metrics, methods)) %>% t)) %>% map(~.[c(1,2,9,10,7,8,3,4,5,6),])
sim1_st_sd <- sim1_st %>% map(~(apply(.,1,sd) %>% matrix(nrow = 6, ncol = 10, dimnames = list(metrics, methods)) %>% t)) %>% map(~.[c(1,2,9,10,7,8,3,4,5,6),])
sim2_sl_sd <- sim2_sl %>% map(~(apply(.[,-c1],1,sd, na.rm = T) %>% matrix(nrow = 6, ncol = 10, dimnames = list(metrics, methods)) %>% t)) %>% map(~.[c(1,2,9,10,7,8,3,4,5,6),])
sim2_st_sd <- sim2_st %>% map(~(apply(.[,-c2],1,sd, na.rm = T) %>% matrix(nrow = 6, ncol = 10, dimnames = list(metrics, methods)) %>% t)) %>% map(~.[c(1,2,9,10,7,8,3,4,5,6),])
t2_sd <- rbind(sim1_sl_sd[[1]], sim1_st_sd[[1]], sim2_sl_sd[[1]], sim2_st_sd[[1]])[, c(2,3,4,1)] %>% round(2)
t2 <- data.frame(method = rownames(t2_mean),
                 tuning = rep(c("AIC", "BIC"), 20),
                 sensitivity = paste0(t2_mean[,"sensitivity"], "(", t2_sd[,"sensitivity"], ")"),
                 specifity = paste0(t2_mean[,"specifity"], "(", t2_sd[,"specifity"], ")"),
                 l2loss = paste0(t2_mean[,"l2loss"], "(", t2_sd[,"l2loss"], ")"),
                 parameter =  paste0(t2_mean[,"Tuning parameter"], "(", t2_sd[,"Tuning parameter"], ")"))
t2
xtable(t2)

# 10var(beta X) = var(epsilon)
c1 <- c(89, 102, 127, 187, 384, 421, 503, 577, 850,944)
c2 <- c(362, 892, 426, 630, 295)
sim1_sl_mean <- sim1_sl %>% map(~(apply(.,1,mean) %>% matrix(nrow = 6, ncol = 10, dimnames = list(metrics, methods)) %>% t))  %>% map(~.[c(1,2,9,10,7,8,3,4,5,6),])
sim1_st_mean <- sim1_st %>% map(~(apply(.,1,mean) %>% matrix(nrow = 6, ncol = 10, dimnames = list(metrics, methods)) %>% t)) %>% map(~.[c(1,2,9,10,7,8,3,4,5,6),])
sim2_sl_mean <- sim2_sl %>% map(~(apply(.[,-c1],1,mean, na.rm = T) %>% matrix(nrow = 6, ncol = 10, dimnames = list(metrics, methods)) %>% t)) %>% map(~.[c(1,2,9,10,7,8,3,4,5,6),])
sim2_st_mean <- sim2_st %>% map(~(apply(.[,-c2],1,mean, na.rm = T) %>% matrix(nrow = 6, ncol = 10, dimnames = list(metrics, methods)) %>% t)) %>% map(~.[c(1,2,9,10,7,8,3,4,5,6),])
t3_mean <- rbind(sim1_sl_mean[[2]], sim1_st_mean[[2]], sim2_sl_mean[[2]], sim2_st_mean[[2]])[, c(2,3,4,1)] %>% round(2)
sim1_sl_sd <- sim1_sl %>% map(~(apply(.,1,sd) %>% matrix(nrow = 6, ncol = 10, dimnames = list(metrics, methods)) %>% t)) %>% map(~.[c(1,2,9,10,7,8,3,4,5,6),])
sim1_st_sd <- sim1_st %>% map(~(apply(.,1,sd) %>% matrix(nrow = 6, ncol = 10, dimnames = list(metrics, methods)) %>% t)) %>% map(~.[c(1,2,9,10,7,8,3,4,5,6),])
sim2_sl_sd <- sim2_sl %>% map(~(apply(.[,-c1],1,sd, na.rm = T) %>% matrix(nrow = 6, ncol = 10, dimnames = list(metrics, methods)) %>% t)) %>% map(~.[c(1,2,9,10,7,8,3,4,5,6),])
sim2_st_sd <- sim2_st %>% map(~(apply(.[,-c2],1,sd, na.rm = T) %>% matrix(nrow = 6, ncol = 10, dimnames = list(metrics, methods)) %>% t)) %>% map(~.[c(1,2,9,10,7,8,3,4,5,6),])
t3_sd <- rbind(sim1_sl_sd[[1]], sim1_st_sd[[2]], sim2_sl_sd[[2]], sim2_st_sd[[2]])[, c(2,3,4,1)] %>% round(2)
t3 <- data.frame(method = rownames(t3_mean),
                 tuning = rep(c("AIC", "BIC"), 20),
                 sensitivity = paste0(t3_mean[,"sensitivity"], "(", t3_sd[,"sensitivity"], ")"),
                 specifity = paste0(t3_mean[,"specifity"], "(", t3_sd[,"specifity"], ")"),
                 l2loss = paste0(t3_mean[,"l2loss"], "(", t3_sd[,"l2loss"], ")"),
                 parameter =  paste0(t3_mean[,"Tuning parameter"], "(", t3_sd[,"Tuning parameter"], ")"))
t3
xtable(t3)

# Scenario 1 ----------
# sim1_sl <- readRDS("sim1_sl.rds")[[1]]
# sensitivity <- map_dfr(1:10, function(k){
#   data.frame(value = sim1_sl[seq(2, 60, by = 6)[k], ],
#              method = rep(c("CT-LASSO", "TASSO", "CLASSO", "TFL-2", "CT-LASSO-p"), each = 2)[k],
#              tuning = rep(c("AIC", "BIC"), 5)[k])
# })
# p1_sl_sen <- ggplot(sensitivity) +
#   geom_boxplot(aes(method, value, color = tuning)) + 
#   theme_bw() +
#   labs(x = "", y = "", color = "Tuning") +
#   coord_flip() +
#   theme(axis.text = element_text(size = 12, face = "bold")) +
#   scale_x_discrete(limits = c("CLASSO", "TASSO", "TFL-2", "CT-LASSO-p",  "CT-LASSO")) +
#   theme(legend.position = "bottom", legend.text = element_text(size = 12), legend.title = element_text(size = 12))
#   
# 
# specificity <- map_dfr(1:10, function(k){
#   data.frame(value = sim1_sl[seq(3, 60, by = 6)[k], ],
#              method = rep(c("CT-LASSO", "TASSO", "CLASSO", "TFL-2", "CT-LASSO-p"), each = 2)[k],
#              tuning = rep(c("AIC", "BIC"), 5)[k])
# })
# p1_sl_spe <- ggplot(specificity) +
#   geom_boxplot(aes(method, value, color = tuning)) + 
#   theme_bw() +
#   labs(x = "", y = "", color = "Tuning") +
#   coord_flip() +
#   theme(axis.text = element_text(size = 12, face = "bold")) +
#   scale_x_discrete(limits = c("CLASSO", "TASSO", "TFL-2", "CT-LASSO-p",  "CT-LASSO")) +
#   theme(legend.position = "bottom", legend.text = element_text(size = 12), legend.title = element_text(size = 12))
# 
# 
# l2loss <- map_dfr(1:10, function(k){
#   data.frame(value = sqrt(sim1_sl[seq(4, 60, by = 6)[k], ]),
#              method = rep(c("CT-LASSO", "TASSO", "CLASSO", "TFL-2", "CT-LASSO-p"), each = 2)[k],
#              tuning = rep(c("AIC", "BIC"), 5)[k])
# })
# p1_sl_l2l <- ggplot(l2loss) +
#   geom_boxplot(aes(method, value, color = tuning)) + 
#   theme_bw() +
#   labs(x = "", y = "", color = "Tuning") +
#   coord_flip() +
#   theme(axis.text = element_text(size = 12, face = "bold")) +
#   scale_x_discrete(limits = c("CLASSO", "TASSO", "TFL-2", "CT-LASSO-p",  "CT-LASSO")) +
#   theme(legend.position = "bottom", legend.text = element_text(size = 12), legend.title = element_text(size = 12))
# 
# tuning_param <- map_dfr(1:4, function(k){
#   data.frame(value = sqrt(sim1_sl[c(1,7,49,55)[k], ]),
#              method = rep(c("CT-LASSO", "CT-LASSO-p"), each = 2)[k],
#              tuning = rep(c("AIC", "BIC"), 5)[k])
# })
# p1_sl_tun <- ggplot(tuning_param) +
#   geom_boxplot(aes(method, value, color = tuning)) + 
#   theme_bw() +
#   labs(x = "", y = "", color = "Tuning") +
#   coord_flip() +
#   theme(axis.text = element_text(size = 12, face = "bold")) +
#   scale_x_discrete(limits = rev) +
#   theme(legend.position = "bottom", legend.text = element_text(size = 12), legend.title = element_text(size = 12))
# 
# p1_sl <- plot_grid(p1_sl_sen, p1_sl_spe, p1_sl_l2l, p1_sl_tun, nrow = 1)
# p1_sl <- plot_grid(ggdraw()+draw_label("Scenario 1", angle = 90, fontface = 'bold', size = 20), p1_sl, rel_widths = c(0.02,1))
# p1_sl <- plot_grid(p1_sl, ggdraw()+draw_line(x=c(0,1), y = c(0.5,0.5)), ncol = 1, rel_heights = c(1, 0.01))
# # save_plot(p1_sl, filename = "p1_sl.png", base_asp = 4)
# 
# # Scenario 2 ---------------
# sim1_st <- readRDS("sim1_st.rds")[[1]]
# sensitivity <- map_dfr(1:10, function(k){
#   data.frame(value = sim1_st[seq(2, 60, by = 6)[k], ],
#              method = rep(c("CT-LASSO", "TASSO", "CLASSO", "TFL-2", "CT-LASSO-p"), each = 2)[k],
#              tuning = rep(c("AIC", "BIC"), 5)[k])
# })
# p1_st_sen <- ggplot(sensitivity) +
#   geom_boxplot(aes(method, value, color = tuning)) + 
#   theme_bw() +
#   labs(x = "", y = "", color = "Tuning") +
#   coord_flip() +
#   theme(axis.text = element_text(size = 12, face = "bold")) +
#   scale_x_discrete(limits = c("CLASSO", "TASSO", "TFL-2", "CT-LASSO-p",  "CT-LASSO")) +
#   theme(legend.position = "bottom", legend.text = element_text(size = 12), legend.title = element_text(size = 12))
# 
# 
# specificity <- map_dfr(1:10, function(k){
#   data.frame(value = sim1_st[seq(3, 60, by = 6)[k], ],
#              method = rep(c("CT-LASSO", "TASSO", "CLASSO", "TFL-2", "CT-LASSO-p"), each = 2)[k],
#              tuning = rep(c("AIC", "BIC"), 5)[k])
# })
# p1_st_spe <- ggplot(specificity) +
#   geom_boxplot(aes(method, value, color = tuning)) + 
#   theme_bw() +
#   labs(x = "", y = "", color = "Tuning") +
#   coord_flip() +
#   theme(axis.text = element_text(size = 12, face = "bold")) +
#   scale_x_discrete(limits = c("CLASSO", "TASSO", "TFL-2", "CT-LASSO-p",  "CT-LASSO")) +
#   theme(legend.position = "bottom", legend.text = element_text(size = 12), legend.title = element_text(size = 12))
# 
# 
# l2loss <- map_dfr(1:10, function(k){
#   data.frame(value = sqrt(sim1_st[seq(4, 60, by = 6)[k], ]),
#              method = rep(c("CT-LASSO", "TASSO", "CLASSO", "TFL-2", "CT-LASSO-p"), each = 2)[k],
#              tuning = rep(c("AIC", "BIC"), 5)[k])
# })
# p1_st_l2l <- ggplot(l2loss) +
#   geom_boxplot(aes(method, value, color = tuning)) + 
#   theme_bw() +
#   labs(x = "", y = "", color = "Tuning") +
#   coord_flip() +
#   theme(axis.text = element_text(size = 12, face = "bold")) +
#   scale_x_discrete(limits = c("CLASSO", "TASSO", "TFL-2", "CT-LASSO-p",  "CT-LASSO")) +
#   theme(legend.position = "bottom", legend.text = element_text(size = 12), legend.title = element_text(size = 12))
# 
# tuning_param <- map_dfr(1:4, function(k){
#   data.frame(value = sqrt(sim1_st[c(1,7,49,55)[k], ]),
#              method = rep(c("CT-LASSO", "CT-LASSO-p"), each = 2)[k],
#              tuning = rep(c("AIC", "BIC"), 5)[k])
# })
# p1_st_tun <- ggplot(tuning_param) +
#   geom_boxplot(aes(method, value, color = tuning)) + 
#   theme_bw() +
#   labs(x = "", y = "", color = "Tuning") +
#   coord_flip() +
#   theme(axis.text = element_text(size = 12, face = "bold")) +
#   scale_x_discrete(limits = rev) +
#   theme(legend.position = "bottom", legend.text = element_text(size = 12), legend.title = element_text(size = 12))
# 
# p1_st <- plot_grid(p1_st_sen, p1_st_spe, p1_st_l2l, p1_st_tun, nrow = 1)
# p1_st <- plot_grid(ggdraw()+draw_label("Scenario 2", angle = 90, fontface = 'bold', size = 20), p1_st, rel_widths = c(0.02,1))
# p1_st <- plot_grid(p1_st, ggdraw()+draw_line(x=c(0,1), y = c(0.5,0.5)), ncol = 1, rel_heights = c(1, 0.01))
# # save_plot(p1_st,filename = "p1_st.png", base_asp = 4)
# 
# # Scenario 3 ----------
# sim2_sl <- readRDS("sim2_sl.rds")[[1]]
# sensitivity <- map_dfr(1:10, function(k){
#   data.frame(value = sim2_sl[seq(2, 60, by = 6)[k], ],
#              method = rep(c("CT-LASSO", "TASSO", "CLASSO", "TFL-2", "CT-LASSO-p"), each = 2)[k],
#              tuning = rep(c("AIC", "BIC"), 5)[k])
# })
# p2_sl_sen <- ggplot(sensitivity) +
#   geom_boxplot(aes(method, value, color = tuning)) + 
#   theme_bw() +
#   labs(x = "", y = "", color = "Tuning") +
#   coord_flip() +
#   theme(axis.text = element_text(size = 12, face = "bold")) +
#   scale_x_discrete(limits = c("CLASSO", "TASSO", "TFL-2", "CT-LASSO-p",  "CT-LASSO")) +
#   theme(legend.position = "bottom", legend.text = element_text(size = 12), legend.title = element_text(size = 12))
# 
# 
# specificity <- map_dfr(1:10, function(k){
#   data.frame(value = sim2_sl[seq(3, 60, by = 6)[k], ],
#              method = rep(c("CT-LASSO", "TASSO", "CLASSO", "TFL-2", "CT-LASSO-p"), each = 2)[k],
#              tuning = rep(c("AIC", "BIC"), 5)[k])
# })
# p2_sl_spe <- ggplot(specificity) +
#   geom_boxplot(aes(method, value, color = tuning)) + 
#   theme_bw() +
#   labs(x = "", y = "", color = "Tuning") +
#   coord_flip() +
#   theme(axis.text = element_text(size = 12, face = "bold")) +
#   scale_x_discrete(limits = c("CLASSO", "TASSO", "TFL-2", "CT-LASSO-p",  "CT-LASSO")) +
#   theme(legend.position = "bottom", legend.text = element_text(size = 12), legend.title = element_text(size = 12))
# 
# 
# l2loss <- map_dfr(1:10, function(k){
#   data.frame(value = sqrt(sim2_sl[seq(4, 60, by = 6)[k], ]),
#              method = rep(c("CT-LASSO", "TASSO", "CLASSO", "TFL-2", "CT-LASSO-p"), each = 2)[k],
#              tuning = rep(c("AIC", "BIC"), 5)[k])
# })
# p2_sl_l2l <- ggplot(l2loss[l2loss$value < 5,]) +
#   geom_boxplot(aes(method, value, color = tuning)) + 
#   theme_bw() +
#   labs(x = "", y = "", color = "Tuning") +
#   coord_flip() +
#   theme(axis.text = element_text(size = 12, face = "bold")) +
#   scale_x_discrete(limits = c("CLASSO", "TASSO", "TFL-2", "CT-LASSO-p",  "CT-LASSO")) +
#   theme(legend.position = "bottom", legend.text = element_text(size = 12), legend.title = element_text(size = 12))
# 
# tuning_param <- map_dfr(1:4, function(k){
#   data.frame(value = sqrt(sim2_sl[c(1,7,49,55)[k], ]),
#              method = rep(c("CT-LASSO", "CT-LASSO-p"), each = 2)[k],
#              tuning = rep(c("AIC", "BIC"), 5)[k])
# })
# p2_sl_tun <- ggplot(tuning_param) +
#   geom_boxplot(aes(method, value, color = tuning)) + 
#   theme_bw() +
#   labs(x = "", y = "", color = "Tuning") +
#   coord_flip() +
#   theme(axis.text = element_text(size = 12, face = "bold")) +
#   scale_x_discrete(limits = rev) +
#   theme(legend.position = "bottom", legend.text = element_text(size = 12), legend.title = element_text(size = 12))
# 
# p2_sl <- plot_grid(p2_sl_sen, p2_sl_spe, p2_sl_l2l, p2_sl_tun, nrow = 1)
# p2_sl <- plot_grid(ggdraw()+draw_label("Scenario 3", angle = 90, fontface = 'bold', size = 20), p2_sl, rel_widths = c(0.02,1))
# p2_sl <- plot_grid(p2_sl, ggdraw()+draw_line(x=c(0,1), y = c(0.5,0.5)), ncol = 1, rel_heights = c(1, 0.01))
# # save_plot(p2_sl,filename = "p2_sl.png", base_asp = 4)
# 
# # Scenario 4 ----------
# sim2_st <- readRDS("sim2_st.rds")[[1]]
# sensitivity <- map_dfr(1:10, function(k){
#   data.frame(value = sim2_st[seq(2, 60, by = 6)[k], ],
#              method = rep(c("CT-LASSO", "TASSO", "CLASSO", "TFL-2", "CT-LASSO-p"), each = 2)[k],
#              tuning = rep(c("AIC", "BIC"), 5)[k])
# })
# p2_st_sen <- ggplot(sensitivity) +
#   geom_boxplot(aes(method, value, color = tuning)) + 
#   theme_bw() +
#   labs(x = "", y = "", color = "Tuning") +
#   coord_flip() +
#   theme(axis.text = element_text(size = 12, face = "bold")) +
#   scale_x_discrete(limits = c("CLASSO", "TASSO", "TFL-2", "CT-LASSO-p",  "CT-LASSO"))+
#   theme(legend.position = "bottom", legend.text = element_text(size = 12), legend.title = element_text(size = 12))
# 
# 
# specificity <- map_dfr(1:10, function(k){
#   data.frame(value = sim2_st[seq(3, 60, by = 6)[k], ],
#              method = rep(c("CT-LASSO", "TASSO", "CLASSO", "TFL-2", "CT-LASSO-p"), each = 2)[k],
#              tuning = rep(c("AIC", "BIC"), 5)[k])
# })
# p2_st_spe <- ggplot(specificity) +
#   geom_boxplot(aes(method, value, color = tuning)) + 
#   theme_bw() +
#   labs(x = "", y = "", color = "Tuning") +
#   coord_flip() +
#   theme(axis.text = element_text(size = 12, face = "bold")) +
#   scale_x_discrete(limits = c("CLASSO", "TASSO", "TFL-2", "CT-LASSO-p",  "CT-LASSO")) +
#   theme(legend.position = "bottom", legend.text = element_text(size = 12), legend.title = element_text(size = 12))
# 
# 
# l2loss <- map_dfr(1:10, function(k){
#   data.frame(value = sqrt(sim2_st[seq(4, 60, by = 6)[k], ]),
#              method = rep(c("CT-LASSO", "TASSO", "CLASSO", "TFL-2", "CT-LASSO-p"), each = 2)[k],
#              tuning = rep(c("AIC", "BIC"), 5)[k])
# })
# p2_st_l2l <- ggplot(l2loss[l2loss$value < 10,]) +
#   geom_boxplot(aes(method, value, color = tuning)) + 
#   theme_bw() +
#   labs(x = "", y = "", color = "Tuning") +
#   coord_flip() +
#   theme(axis.text = element_text(size = 12, face = "bold")) +
#   scale_x_discrete(limits = c("CLASSO", "TASSO", "TFL-2", "CT-LASSO-p",  "CT-LASSO")) +
#   theme(legend.position = "bottom", legend.text = element_text(size = 12), legend.title = element_text(size = 12))
# 
# tuning_param <- map_dfr(1:4, function(k){
#   data.frame(value = sqrt(sim2_st[c(1,7,49,55)[k], ]),
#              method = rep(c("CT-LASSO", "CT-LASSO-p"), each = 2)[k],
#              tuning = rep(c("AIC", "BIC"), 5)[k])
# })
# p2_st_tun <- ggplot(tuning_param) +
#   geom_boxplot(aes(method, value, color = tuning)) + 
#   theme_bw() +
#   labs(x = "", y = "", color = "Tuning") +
#   coord_flip() +
#   theme(axis.text = element_text(size = 12, face = "bold")) +
#   scale_x_discrete(limits = rev) +
#   theme(legend.position = "bottom", legend.text = element_text(size = 12), legend.title = element_text(size = 12))
# 
# p2_st <- plot_grid(p2_st_sen, p2_st_spe, p2_st_l2l, p2_st_tun, nrow = 1)
# p2_st <- plot_grid(ggdraw()+draw_label("Scenario 4", angle = 90, fontface = 'bold', size = 20), p2_st, rel_widths = c(0.02,1))
# p2_st <- plot_grid(p2_st, ggdraw()+draw_line(x=c(0,1), y = c(0.5,0.5)), ncol = 1, rel_heights = c(1, 0.01))
# # save_plot(p2_st,filename = "p2_st.png", base_asp = 4)


# header -------
# header <- plot_grid(ggdraw()+draw_label("", fontface = 'bold', size = 20),
#                     ggdraw()+draw_label("Sensitivity", fontface = 'bold', size = 20, hjust = 0.25), 
#                     ggdraw()+draw_label("Specificity", fontface = 'bold', size = 20, hjust = 0.3), 
#                     ggdraw()+draw_label(expression(bold(paste("L"[2], "-loss"))), fontface = 'bold', size = 20), 
#                     ggdraw()+draw_label(expression(bold(paste("Paramter ", eta))), fontface = 'bold', size = 20, hjust = 0.3),
#                     rel_widths = c(0.02,1,1,1,1),nrow = 1)
# header <- plot_grid(header, ggdraw()+draw_line(x=c(0,1), y = c(0.5,0.5)), ncol = 1, rel_heights = c(1, 0.01))
# 
# p <- plot_grid(header, p1_sl, p1_st, p2_sl, p2_st, ncol = 1, rel_heights = c(0.1, 1,1,1,1))
# save_plot(p,filename = "p.png", base_width = 15, base_height = 13)
