setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list=ls())
invisible(sapply(list.files("Funcoes auxiliares",pattern="*.R$",full.names=TRUE, 
                            ignore.case=TRUE),source,.GlobalEnv))
source("Geracao de dados/dst.R")
source("Geracao de dados/mixdst.R")
load("Geracao de dados/n1000type2.Rdata")
library(foreach)
library(doParallel)

# read.table("Geracao de dados/tonedata.txt",h=T, dec = ",", sep = ";") #
dados = as.matrix(data[,-1]) # primeira coluna contem grupos reais

R = 10 # numero de replicacoes
file = "n100type2G2"

cl <- makeCluster(4) 
registerDoParallel(cl)
clusterSetRNGStream(cl, 1)

result = foreach(re = 1:R, .packages = "compiler",.verbose=T) %dopar% {
 enableJIT(3)
 waic_file(file, re, dados, relabel = F)
};stopCluster(cl)

tab = matrix(NA, R, 2)

for(i in 1:R){
 tab[i,] = as.numeric(result[[i]])
}

write.table(tab,"temp.txt", sep = ";", dec = ",", col.names = F, row.names = F);beepr::beep()

