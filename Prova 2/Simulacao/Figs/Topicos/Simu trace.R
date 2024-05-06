setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list=ls())
wd = getwd()
local = paste0(stringr::str_sub(wd,end=-13))
setwd(local)

R = 10 # numero de replicacoes

n = 1000 # tamanho amostral
g = c("G2","","G4") # numero de componentes
type = 1:3 # tipos

windows(500,400)
par(mfrow=c(3,3), mar = c(2.1, 4.5, 3.5, 1.1))
for(j in type){
 for(i in g){
  xx = vector(mode="list",length = R)
  file = paste0("n",n,"type",j,i)
  for(re in 1:R){
   nu.samp = read.table(paste0("Outputs/Paralelo/",file,"/",re, "/param.txt"), h=T)[,1]
   xx[[re]] = coda::mcmc(log(nu.samp))
  }
  xx = coda::mcmc.list(xx)
  coda::traceplot(xx, ylab = expression(log(nu^(q))))
 }
}

