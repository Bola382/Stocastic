# ==============================================================================
# Replicacoes em paralelo do GS utlizando o "Funcao Gibbs sampling.R"
# ==============================================================================

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list=ls())
invisible(sapply(list.files("Funcoes auxiliares",pattern="*.R$",full.names=TRUE, 
                            ignore.case=TRUE),source,.GlobalEnv))
source("Geracao de dados/dst.R")
load("Geracao de dados/dados.Rdata")
source("Gibbs sampling file.R")
library(foreach)
library(doParallel)

dados = as.matrix(data[,-1]) # primeira coluna contem grupos reais

G = 3 # numero de comps
R = 15 # numero de replicacoes
Q =  260000 # iteracoes por replica
burn = 10000 # burnin
thin = 50 # salto
nburn = ceiling((Q-burn)/thin) # numero de amostras pos burn

file = "test2"
dir.create(paste0("Outputs/Paralelo/",file))
sapply(1:R, function(a) dir.create(paste0("Outputs/Paralelo/",file,"/",a),showWarnings = F))
cl <- makeForkCluster(15) # linux
registerDoParallel(cl)
clusterSetRNGStream(cl, 7481)

tempo = foreach(re = 1:R, .packages = "compiler",.verbose=T) %dopar% {
 gibbs(file,re,dados,G,phi=.3,Q,burn,thin)
};stopCluster(cl)
