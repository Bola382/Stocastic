# ==============================================================================
# Gerando amostras da posteriori utilizando Gibbs
# ==============================================================================

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list=ls())

source("Geracao de dados/dst.R")
source("fullnu.R")
source("rtruncnorm.R")
load("Geracao de dados/dados.Rdata")

head(data)

n = nrow(data)
p = ncol(data)-1 # -1 pois a primeira coluna contem os grupos
#obs: p = k+1

y = data[,2] # resp
X = as.matrix(cbind(1,data[,-(1:2)])) # covariaveis com intercepto

# ------------------------
# Inicio do algoritmo
# ------------------------

Q = 10000 # numero de iteracoes de Gibbs
R = 2 # numero de cadeias
phi = 2.8 # desvio padrao do passeio aleatorio de nu
cont = 0 # contador de aceites de MH

# ~~~~~~~~~~~~~~~~
# valores iniciais
# ~~~~~~~~~~~~~~~~

G = 2 # n componentes

# propostas geradas a partir das prioris
set.seed(2)

# para guardar amostras
prob.samp = tau2.samp = Delta.samp = matrix(NA, nrow = Q, ncol = G)
u.samp = t.samp = z.samp = matrix(NA, nrow = Q, ncol = n)
beta.samp = array(NA, dim = c(Q,p,G), dimnames = list(1:Q, 1:p, 1:G))

# ~~~~~~~~~~~~~~~~
# Valores iniciais
# ~~~~~~~~~~~~~~~~

# hiperparametros
xi = rep(1,G) # prob
c = 10 # da beta (c <- sqrt(c) das minhas contas)
eta = 0; omega = 10 # da Delta
r = s = .01 # da tau2

prob.samp[1,] = gtools::rdirichlet(1, alpha = xi) # pesos
beta.samp[1,,] =  matrix(rnorm(p*G,sd=c),ncol=2) # coef reg
tau2.samp[1,] = 1/rgamma(G,r,s)# escala
Delta.samp[1,] = rnorm(G,eta,omega)# forma
alpha.samp = runif(1,.02,.5)# hiperparametro da priori de nu
nu.samp = rexp(1,rate=alpha.samp)# gl

u.samp[1,] = rgamma(n, shape = nu.samp/2, rate = nu.samp/2) # t e u da representacao aumentada
t.samp[1,] = abs(rnorm(n))/sqrt(u.samp[1,])
z.samp[1,] = sample(1:G,n,prob=prob.samp[1,], replace = T) # latente que indica os grupos

# ~~~~~~~~~~~~~~~~~~~
# Barra de progresso
# ~~~~~~~~~~~~~~~~~~~
library(progress)
format = "(:spin) [:bar] :percent [Decorrido: :elapsedfull || Estimado: :eta] Taxa de aceitacao :taxa"
pb = progress_bar$new(format, clear = FALSE, total = Q, complete = "=", incomplete = "-", width = 100)

library(compiler)
enableJIT(3)

set.seed(1)
for(i in 2:Q){
 b = -sqrt(nu.samp[i-1]/pi)*(gamma((nu.samp[i-1]-1)/2)/gamma(nu.samp[i-1]/2))
 U = diag(u.samp[i-1,])
 m = NULL
 #m = table(z.samp[i-1,]) # n de amostras por comp
 
 # ----------------------
 # Atualizando alpha
 # ----------------------
 
 # gama truncada em (0.02,0.5)
 alpha.samp[i] = Runuran::urgamma(1, shape = 2, scale = 1/nu.samp[i-1], lb = .02, ub = .5)
 
 #print(u.samp[i-1,1:10])
 
 # atualizando pametros por componente
 for(j in 1:G){
  index = which(z.samp[i-1,]==j)
  mj = length(index)
  m[j] = mj
  Xj = X[index,] # dim = mjXp
  Uj = U[index,index] # dim = mjXmj
  yj = y[index] # dim = mjX1
  tj = t.samp[i-1, index] # dim = mjX1
  
  if(mj==1){
   XU = as.matrix(Xj)%*%Uj # t(Xj)%*%Uj
  }else{
   XU = crossprod(Xj,Uj) # t(Xj)%*%Uj
  } 
  aux_y = yj-(b+tj)*Delta.samp[i-1,j] 
  sigma_inv = XU%*%Xj + diag(tau2.samp[i-1,j]/c^2,nrow = p)
  
  # ----------------------
  # Atualizando beta
  # ----------------------
  #s_sigma = solve(sigma_inv)
  mmm = solve(sigma_inv,XU%*%aux_y) #crossprod(s_sigma,XU%*%aux_y)
  ppp = tau2.samp[i-1,j]*solve(sigma_inv)#*s_sigma
  # [iter,coefs,componente]
  beta.samp[i,,j] = MASS::mvrnorm(1, mu = mmm, 
                                     Sigma = ppp)
  #if(j==1){cat(i,"\t", beta.samp[i,,j], "\t med",mmm[1],"\t var", ppp[1,1], "\n")}
  # ----------------------
  # Atualizando tau2
  # ----------------------
  S3 = sum(u.samp[i-1,index]*(yj-Xj%*%beta.samp[i,,j]-(b+tj)*Delta.samp[i-1,j])^2)
  
  tau2.samp[i,j] = LearnBayes::rigamma(1, mj/2+r, (S3+s)/2)
  #MCMCpack::rinvgamma(1, mj/2+r, 2/(S3+s))
  #1/rgamma(1,shape = m[j]/2+r, rate = (S3 + s)/2)
 
  # ----------------------
  # Atualizando Delta
  # ----------------------
  S1 = sum(u.samp[i-1,index]*(yj-Xj%*%beta.samp[i,,j])*(b+tj))
  S2 = sum(u.samp[i-1,index]*(b+tj)^2)
  denom = omega^2*S2 + tau2.samp[i,j]
  
  # aqui tomamos eta = 0 para evitar mais conta, lembrar de 
  # alterar a media caso eta != 0.
  Delta.samp[i,j] = rnorm(1, omega^2*S1/denom, omega*sqrt(tau2.samp[i,j]/denom))
 }
 
 # ----------------------
 # Atualizando prob
 # ----------------------
 
 prob.samp[i,] = gtools::rdirichlet(1, alpha = xi + m)
 
 # ------------------------------
 # atualizando variaveis latentes
 # ------------------------------
 aux_mu = X%*%beta.samp[i,,]+matrix(rep(b*Delta.samp[i,],n), nrow = n, ncol = G, byrow = T)
 aux_sig = sqrt(tau2.samp[i,]+Delta.samp[i,]^2)
 aux_lam = Delta.samp[i,]/sqrt(tau2.samp[i,])
 
 for(k in 1:n){
  # ----------------------
  # Atualizando Z
  # ----------------------
  aux_prob = prob.samp[i,]*sapply(1:G,function(a) dst(y[k], aux_mu[k,a], aux_sig[a], aux_lam[a], nu.samp[i-1]))
  z.samp[i,k] = sample(1:G,1, prob = aux_prob) # prob eh normalizado internamente
  
  # ----------------------
  # Atualizando T
  # ----------------------
  aux_yxbeta = y[k]-aux_mu[k,z.samp[i,k]]
  aux_denom = Delta.samp[i,z.samp[i,k]]^2 + tau2.samp[i,z.samp[i,k]]
  aux_mean = (aux_yxbeta*Delta.samp[i,z.samp[i,k]])/aux_denom
  aux_sd = sqrt(tau2.samp[i,z.samp[i,k]]/(u.samp[i-1,k]*aux_denom))
  # media antes do truncamento muito distante do lb
  
  aux_truncnorm = rtruncnorm(1,-aux_mean/aux_sd)$amostra
  
  t.samp[i,k] = aux_mean+aux_sd*aux_truncnorm
  # ----------------------
  # Atualizando U
  # ----------------------
  D1 =(aux_yxbeta-t.samp[i,k]*Delta.samp[i,z.samp[i,k]])^2/(2*tau2.samp[i,z.samp[i,k]]) 
  u.samp[i,k] = rgamma(1, shape = nu.samp[i-1]/2+1, rate = D1 + (t.samp[i,k]^2+nu.samp[i-1])/2)
 }
 # ----------------------
 # Atualizando nu
 # ----------------------
 
 # Passo de MH
 prop = rlnorm(1, log(nu.samp[i-1]), phi)
 
 aceit = min(1,exp(fullnu(prop,alpha.samp[i],prob.samp[i,],y,aux_mu,aux_sig,aux_lam) + log(prop) - fullnu(nu.samp[i-1],alpha.samp[i],prob.samp[i,],y,aux_mu,aux_sig,aux_lam) - log(nu.samp[i-1])))
 if(aceit > 1 | aceit <0){stop("Problema no passo MH")}
 if(runif(1)<=aceit){
  nu.samp[i] = prop
  cont = cont + 1
 }else{
  nu.samp[i] = nu.samp[i-1]
 }
 pb$tick(tokens = list(taxa = paste0(formatC(cont/i * 100,2,format="f"),"%")))
};beepr::beep()

# aceitacao do MH deve estar entre 0.234 e 0.44

library(coda)

# [iter,coefs,componente]
traceplot(mcmc(beta.samp[,1,1]))
traceplot(mcmc(beta.samp[,,2]))
acf(beta.samp[,,1])
acf(beta.samp[,,2])
# 
# traceplot(mcmc(tau2.samp[200:Q,2]))
# traceplot(mcmc(tau2.samp))
# acf(beta.samp[,,1])
# acf(beta.samp[,,2])
# 
# traceplot(mcmc(nu.samp[200:Q]))
# traceplot(mcmc(nu.samp))
# acf(nu.samp)
# 
# i = which.max(beta.samp[,1,1])
# j = 1
# 
# b = -sqrt(nu.samp[i-1]/pi)*gamma((nu.samp[i-1]-1)/2)/gamma(nu.samp[i-1]/2)
# 
# 
# U = diag(u.samp[i-1,])
# m = NULL
# index = which(z.samp[i-1,]==j)
# mj = length(index)
# m[j] = mj
# Xj = X[index,] # dim = mjXp
# Uj = U[index,index] # dim = mjXmj
# yj = y[index] # dim = mjX1
# tj = t.samp[i-1, index] # dim = mjX1
# 
# if(mj==1){
#  XU = as.matrix(Xj)%*%Uj # t(Xj)%*%Uj
# }else{
#  XU = crossprod(Xj,Uj) # t(Xj)%*%Uj
# }
# aux_y = yj-(b+tj)*Delta.samp[i-1,j]
# sigma_inv = XU%*%Xj + diag(tau2.samp[i-1,j]/c^2,nrow = p)
# 
# solve(sigma_inv)%*%sigma_inv
# 
# lu_decomp = matrixcalc::lu.decomposition(sigma_inv)
# solve(lu_decomp$U, solve(lu_decomp$L, diag(1, nrow = dim(lu_decomp$L))))%*%sigma_inv
# 
# solve(sigma_inv,crossprod(t(XU),aux_y))
# tau2.samp[i-1,j]*
# 
# MASS::mvrnorm(1, mu = solve(sigma_inv,XU%*%aux_y),
#                                  Sigma = tau2.samp[i-1,j]*solve(sigma_inv))
# k=1
# aux_yxbeta = y[k]-aux_mu[k,z.samp[i,k]]
# aux_denom = Delta.samp[i,z.samp[i,k]]^2 + tau2.samp[i,z.samp[i,k]]
# aux_mean = (aux_yxbeta*Delta.samp[i,z.samp[i,k]])/aux_denom
# aux_sd = sqrt(tau2.samp[i,z.samp[i,k]]/(u.samp[i-1,k]*aux_denom))
# -aux_mean/aux_sd
