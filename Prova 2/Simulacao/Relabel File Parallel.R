# ==============================================================================
# Aplica "Re-labeling File.R" em paralelo
# ==============================================================================

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list=ls())
source("Re-labeling File.R")
library(foreach)
library(doParallel)

R = 10 # numero de replicacoes
file = "n1000type2G2"

cl <- makeCluster(4) 
registerDoParallel(cl)
clusterSetRNGStream(cl, 1)

result = foreach(re = 1:R, .packages = "compiler",.verbose=T) %dopar% {
 relabel(file,re, try = 100)
};stopCluster(cl);beepr::beep();result

# write.table(unlist(result),"temp.txt", sep = ";", dec = ",", col.names = F, row.names = F)