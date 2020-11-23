library(data.tree)
setwd("~/Dropbox/research/graphical-model/hierarchical-lasso/compositional-hierarchical-tree-regression/")
my.roi <- readRDS("data-analysis/tree-structure.rds")
q <- ncol(my.roi)
my.roi$pathString <- paste("Whole brain", my.roi$`level1-type1`, sep = "/")
for(j in 1:nrow(my.roi)){
  string <- my.roi$pathString[j]
  for(k in (q-3):1){
    if(my.roi[j,k] != my.roi[j,k+1]){
      string <- paste(string, my.roi[j,k], sep = "/")
    }
  } 
  my.roi$pathString[j] <- gsub('.lvl[0-9]', '', string)
}

tree.structure <- as.Node(my.roi)
print(tree.structure)
plot(tree.structure)
library(igraph)
plot(as.igraph(tree.structure, directed = TRUE, direction = "climb"))
colorVector <- c(rep("red", 10), rep("black", 10))
jsarray <- paste0('["', paste(colorVector, collapse = '", "'), '"]')
nodeStrokeJS <- JS(paste0('function(d, i) { return ', jsarray, '[i]; }'))
radialNetwork(ToListExplicit(tree.structure, unname = TRUE), nodeStroke = nodeStrokeJS)


library(networkD3)
diagonalNetwork(ToListExplicit(tree.structure, unname = TRUE), fontSize = 10, opacity = 0.9, width = 800, height=2000)
radialNetwork(ToListExplicit(tree.structure, unname = TRUE), fontSize = 9, opacity = 0.9, width = 800, height=800)


CTASSO_BIC <- readRDS("data-analysis/CTASSO-BIC.rds")
CTASSO_AIC <- readRDS("data-analysis/CTASSO-AIC.rds")
CTASSO_CV <- readRDS("data-analysis/CTASSO-CV.rds")
