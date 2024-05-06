setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list=ls())
wd = getwd()
local = paste0(stringr::str_sub(wd,end=-13))
setwd(local)

R = 10 # numero de replicacoes

n = 1000 # tamanho amostral
g = c("G2","","G4") # numero de componentes
type = 1:3 # tipos

for(j in type){
 for(i in g){
  file = paste0("n",n,"type",j,i)
  for(re in 1:R){
   caminho = paste0("Outputs/Paralelo/",file,"/",re,"/")
   
   beta.samp = read.table(paste0(caminho,"beta.txt"), h =T)
   tau2.samp = read.table(paste0(caminho,"tau2.txt"), h =T)
   Delta.samp = read.table(paste0(caminho,"Delta.txt"), h =T)
   prob.samp = read.table(paste0(caminho,"prob.txt"), h =T)
   z.samp = read.table(paste0(caminho,"z.txt"), h =T)
   ut = read.table(paste0(caminho,"ut.txt"), h =T)
   param = read.table(paste0(caminho,"param.txt"), h =T)
   
   Q = nrow(beta.samp)
   write.table(beta.samp[100:Q,],paste0(caminho,"beta.txt"),row.names = F)
   write.table(tau2.samp[100:Q,],paste0(caminho,"tau2.txt"),row.names = F)
   write.table(Delta.samp[100:Q,],paste0(caminho,"Delta.txt"),row.names = F)
   write.table(prob.samp[100:Q,],paste0(caminho,"prob.txt"),row.names = F)
   write.table(z.samp[100:Q,],paste0(caminho,"z.txt"),row.names = F)
   write.table(ut[100:Q,],paste0(caminho,"ut.txt"),row.names = F)
   write.table(param[100:Q,],paste0(caminho,"param.txt"),row.names = F)
  }
 }
}
