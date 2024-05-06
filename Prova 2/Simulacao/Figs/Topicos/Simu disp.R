setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list=ls())
wd = getwd()
local = paste0(stringr::str_sub(wd,end=-13),"Geracao de dados")
setwd(local)

n = c(100,1000) # tamanhos amostrais
type = 1:3 # tipos

windows(300,350)
par(mfrow=c(3,2), mar = c(1.1, 4, 4.1, 1.1))
for(j in type){
 for(i in n){
  load(paste0("n",i,"type",j,".RData"))
  if(j == 3 & (i == 100 || i == 1000)){par(mar = c(4.1, 4, 1.5, 1.1))}
  col = sapply(data$comp, function(a) if(a==1){rgb(1,0,0)}else{if(a==2){rgb(0,1,0)}else{rgb(0,0,1)}})
  plot(data$x1, data$resp, col = col, cex = .9, pch = 16, 
       xlab = ifelse(j==1|j==2,"","x"), ylab = ifelse(i==1000,"","y"), 
       xlim = c(-3.5,3.5), ylim = c(min(data$resp)-5,max(data$resp)+ifelse(i==1000 & j == 1,100,5)))
  apply(beta,2,function(a) abline(a, lty=2)) 
  if(i==1000 & j == 1){
   legend("topright", legend = paste("Componente", 1:3), pch = 16, 
          col = c(rgb(1,0,0),
                  rgb(0,1,0),
                  rgb(0,0,1)))
  }
  rm(list =  setdiff(ls(), c("i","j","n","type")))
 }
}

windows(350,300)
par(mfrow=c(2,3), mar = c(1.1, 4, 4.1, 1.1))

for(i in n){
  for(j in type){
   load(paste0("n",i,"type",j,".RData"))
   if(i == 100 & (j == 1 || j == 2 || j==3)){par(mar = c(4.1, 4, 1.5, 1.1))}
   col = sapply(data$comp, function(a) if(a==1){rgb(1,0,0)}else{if(a==2){rgb(0,1,0)}else{rgb(0,0,1)}})
   plot(data$x1, data$resp, col = col, cex = .9, pch = 16, 
        xlab = ifelse(i==100,"","x"), ylab = ifelse(j==2 || j==3,"","y"), 
        xlim = c(-3.5,3.5), ylim = c(min(data$resp)-5,max(data$resp)+ifelse(i==100 & j == 3,20,5)))
   apply(beta,2,function(a) abline(a, lty=2)) 
   if(i==100 & j == 3){
    legend("topright", legend = paste("Componente", 1:3), pch = 16, 
           col = c(rgb(1,0,0),
                   rgb(0,1,0),
                   rgb(0,0,1)))
  }
  rm(list =  setdiff(ls(), c("i","j","n","type")))
  }
}
