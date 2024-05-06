setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list=ls())
invisible(sapply(list.files("Funcoes auxiliares",pattern="*.R$",full.names=TRUE, 
                            ignore.case=TRUE),source,.GlobalEnv))
source("Geracao de dados/dst.R")
source("Geracao de dados/mixpst.R")
load("Geracao de dados/n100type3.Rdata")
library(foreach)
library(doParallel)

# read.table("Geracao de dados/tonedata.txt",h=T, dec = ",", sep = ";") 
dados = as.matrix(data[,-1]) # primeira coluna contem grupos reais

R = 10 # numero de replicacoes
file = "n100type3G4"

cl <- makeCluster(4) 
registerDoParallel(cl)
clusterSetRNGStream(cl, 1)

result = foreach(re = 1:R, .packages = c("compiler","sn"),.verbose=T) %dopar% {
 enableJIT(3)
 qqplot_file(file, re, dados)
};stopCluster(cl)

for(i in 1:R){
 res = result[[i]]
 write.table(res,paste0("Outputs/Paralelo/",file,"/",i,"/","res.txt"), row.names = F)
};beepr::beep()

