setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list=ls())

file = "RealG3"
R = 10
tab = matrix(NA,R,2)

for(i in 1:R){
 tab[i,] =  as.numeric(read.table(paste0("Outputs/Paralelo/",file,"/",i, "/monitor.txt"), h=T))
}

write.table(tab,"temp.txt", sep = ";", dec = ",", col.names = F, row.names = F)
