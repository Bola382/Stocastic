setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list=ls())

R = 10 # numero de replicacoes
file = "n1000type3G2"

xx = vector(mode="list",length = R)

for(re in 1:R){
 xx[[re]] = coda::mcmc(log(read.table(paste0("Outputs/Paralelo/",file,"/",re, "/param.txt"), h=T)[,1]))
}

xx = coda::mcmc.list(xx)
coda::traceplot(xx, xlab = "Iterações", ylab = expression(log(nu^(q))))
coda::gelman.diag(xx)
