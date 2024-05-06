setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list=ls())

R = 10 # numero de replicacoes

g = c("G2","","G4") # numero de componentes

type = 3 # tipos
n = 100 # tamanho amostral

file = paste0("n",n,"type",type,g[2])
file2 = paste0("n",n,"type",type)

#file = paste0("Real",g[1])

load(paste0("Geracao de dados/",file2,".RData"))
dados =  data[,-1] # primeira coluna tem as alocacoes
#dados = read.table("Geracao de dados/tonedata.txt",h=T, dec = ",", sep = ";")

moda = function(vec){
 unique(sort(vec))[which.max(table(vec))]
}

contagem = function(vec,categ){
 sapply(1:categ, function(a) sum(vec==a))
}

z_list = matrix(NA, nrow = R, ncol = nrow(dados)) # para classificacao por replica

mcmc.samp = NULL # tudo junto
z.mcmc = NULL  # classificiacao tudo junto
compiler::enableJIT(3)
for(repl in 1:R){
 # recuperando parametros por rep
 caminho = paste0("Outputs/Paralelo/",file,"/",repl,"/corrigido/")
 beta.samp = read.table(paste0(caminho,"beta.txt"), h =T)
 tau2.samp = read.table(paste0(caminho,"tau2.txt"), h =T)
 Delta.samp = read.table(paste0(caminho,"Delta.txt"), h =T)
 prob.samp = read.table(paste0(caminho,"prob.txt"), h =T)
 param = read.table(paste0(caminho,"param.txt"), h =T)
 nu.samp = param$nu;rm("param")
 
 lambda.samp = Delta.samp/sqrt(tau2.samp)
 sigma2.samp = tau2.samp+Delta.samp^2
 z.samp = read.table(paste0(caminho,"z.txt"), h =T)
 
 Q = nrow(beta.samp)
 #z_tmp = apply(z.samp,2,moda) # classificacao pelo maximo a posteriori
 beta_tmp = matrix(as.matrix(beta.samp[Q,]),nrow = ncol(dados), ncol = ncol(beta.samp)/ncol(dados))
 
 # rerotulando baseado nas retas de reg
 ordem = crossprod(c(1,max(as.numeric(dados[,2]))),beta_tmp)
 novaordem = order(ordem, decreasing = T)
 
 #z_list[repl,] = sapply(z_tmp, function(a) which(novaordem==a)) # reordena
 z_ok = t(apply(z.samp,1, function(a) sapply(a, function(b) which(novaordem==b))))
 
 
 G = ncol(tau2.samp)
 prob_ok = prob.samp[,novaordem]; colnames(prob_ok) = paste0("p_",1:G)
 sigma2_ok = sigma2.samp[,novaordem]; colnames(sigma2_ok) = paste0("sigma2_",1:G)
 lambda_ok = lambda.samp[,novaordem]; colnames(lambda_ok) = paste0("lambda_",1:G)
 
 
 # organizando beta.samp em um array
 n = ncol(z.samp)
 p = ncol(beta.samp)/G # numero de betas (covs + 1)
 beta.aux = array(NA, dim = c(Q,p,G), dimnames = list(1:Q,1:p,1:G))
 beta.samp = as.matrix(beta.samp)
 
 for(i in 1:Q){
  qual = c(1:(p-1),0)
  for(j in 1:p){
   quais = which(1:ncol(beta.samp) %% p == qual[j])
   beta.aux[i,j,] = beta.samp[i,quais]
  }
 }
 beta.aux = beta.aux[,,novaordem]
 # convertendo beta_ok em uma matriz Q x (p*G)
 beta_ok = matrix(NA, nrow = Q, ncol = p*G, dimnames = list(1:Q,paste0(rep(paste0("beta",1:p),G),"_",rep(1:G,each=p))))
 
 for(ii in 1:Q){
  beta_ok[ii,] = c(beta.aux[ii,,])
 }
 
 mcmc.samp = rbind(mcmc.samp,cbind(prob_ok,beta_ok,sigma2_ok,lambda_ok,nu.samp))
 z.mcmc = rbind(z.mcmc,z_ok)
}
#colnames(mcmc.samp) = c(paste0("p_",1:G),paste0(rep(paste0("beta",1:p),G),"_",rep(1:G,each=p)),
#                        paste0("sigma2_",1:G),paste0("lambda_",1:G),"nu")
colMeans(mcmc.samp)
#apply(z.mcmc,2,moda)

# ----------------------------
# por replica
# ----------------------------

# assume que o numero de G nao ultrapassa 3
# col = apply(z_list, 2,function(a) contagem(a,3)/R)
# col = apply(col,2,function(a) rgb(a[1],a[2],a[3]))
# 
# # plot
# plot(dados[,2],dados[,1], col=col, pch = 16, cex = .9, 
#      ylab = "y", xlab = "x")

# ----------------------------
# tudo junto
# ----------------------------

# assume que o numero de G nao ultrapassa 3
Q = nrow(z.mcmc)
beta.hat = colMeans(mcmc.samp)[(G+1):(G+p*G)]
beta.hat = matrix(beta.hat,nrow = p, ncol = G)

col = apply(z.mcmc, 2,function(a) contagem(a,3)/Q)
col = apply(col,2,function(a) rgb(a[1],a[2],a[3]))

# plot
# windows(350,300)
# par(mfrow=c(2,3), mar = c(1.1, 4, 4.1, 1.1))
plot(dados[,2],dados[,1], col=col, pch = 16, cex = .9,
     ylab = "y", xlab = "x")
apply(beta.hat,2,function(a) abline(a))
apply(beta,2,function(a) abline(a, lty=2))
