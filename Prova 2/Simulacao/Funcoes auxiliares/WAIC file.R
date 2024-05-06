# ==================================================================
# Calculo do WAIC paralelizado
# ===================================================================
# file: nome da pasta em Outputs/Paralelo onde estao salvos os resultados
# repl: qual replicacao escolher
# dados: matriz de dados, coluna 1 respostas e restantes covariaveis

waic_file = function(file,repl,dados,relabel=T){
 caminho = paste0("Outputs/Paralelo/",file,"/",repl,ifelse(relabel,"/corrigido/","/"))

 beta.samp = read.table(paste0(caminho,"beta.txt"), h =T)
 tau2.samp = read.table(paste0(caminho,"tau2.txt"), h =T)
 Delta.samp = read.table(paste0(caminho,"Delta.txt"), h =T)
 prob.samp = read.table(paste0(caminho,"prob.txt"), h =T)
 param = read.table(paste0(caminho,"param.txt"), h =T)
 nu.samp = param$nu;rm("param")
 
 lambda.samp = Delta.samp/sqrt(tau2.samp)
 sigma.samp = sqrt(tau2.samp+Delta.samp^2)
 
 waic(dados, beta.samp, prob.samp, sigma.samp, lambda.samp, nu.samp)
}