setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list=ls())
wd = getwd()
local = paste0(stringr::str_sub(wd,end=-13))
setwd(local)

R = 10 # numero de replicacoes

n = 1000 # tamanho amostral
g = c("G2","","G4") # numero de componentes
type = 1:3 # tipos

quantis = qnorm((1:n-.5)/n)

windows(500,400)
par(mfrow=c(3,3), mar = c(2.1, 4.5, 3.5, 1.1))
for(j in type){
 for(i in g){
  if(j==1 & i=="G2"){
   plot(1,1,type="n", xlab = "Quantis observados", ylab = "Quantis teóricos")
  }else{
   res = vector(mode="list", length = R)
   file = paste0("n",n,"type",j,i)
   for(re in 1:R){
    res[[re]] = read.table(paste0("Outputs/Paralelo/",file,"/",re, "/res.txt"), h=T)
   }
   tmp.min = lapply(res, function(a) apply(a,2,min))
   tmp.max = lapply(res, function(a) apply(a,2,max))
   xmin = min(sapply(1:R, function(a) tmp.min[[a]][2]))
   xmax = max(sapply(1:R, function(a) tmp.max[[a]][3]))
   
   plot(1,1,type="n", xlab = "Quantis observados", ylab = "Quantis teóricos",
        ylim = c(-3,3), xlim = c(xmin-1.5,xmax+1.5))
   set.seed(1)
   invisible(sapply(1:R, function(a) points(res[[a]][,1],quantis, 
                                            pch = 16, cex = .9,col = a)))
   abline(a=0,b=1)
   invisible(sapply(1:R, function(a) lines(res[[a]][,2],quantis, 
                                           lty = 2,col = a)))
   invisible(sapply(1:R, function(a) lines(res[[a]][,3],quantis, 
                                           lty = 2,col = a)))
   }
 }
}

