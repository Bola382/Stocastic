# ==============================================================================
# Gerando amostras da posteriori utilizando Telescoping
# ==============================================================================

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list=ls())
sapply(list.files("Funcoes auxiliares",pattern="*.R$",full.names=TRUE, 
                  ignore.case=TRUE),source,.GlobalEnv)
source("Geracao de dados/dst.R")
source("Funcoes auxiliares/fullnu.R")
source("Funcoes auxiliares/rtruncnorm.R")
load("Geracao de dados/dados.Rdata")

head(data)

n = nrow(data)
p = ncol(data)-1 # -1 pois a primeira coluna contem os grupos
#obs: p = k+1

y = data[,2] # resp
X = as.matrix(cbind(1,data[,-(1:2)])) # covariaveis com intercepto

Q = 50000 # numero de iteracoes

# ~~~~~~~~~~~~~~~~
# valores iniciais
# ~~~~~~~~~~~~~~~~

G = 3 # n componentes

# propostas geradas a partir das prioris
set.seed(2)

# para guardar amostras
prob.samp = tau2.samp = Delta.samp = matrix(NA, nrow = Q, ncol = G)
u.samp = t.samp = z.samp = matrix(NA, nrow = Q, ncol = n)
beta.samp = array(NA, dim = c(Q,p,G), dimnames = list(1:Q, 1:p, 1:G))

# hiperparametros
xi = rep(1,G) # prob
c = 10 # da beta (c <- sqrt(c) das minhas contas)
eta = 0; omega = 10 # da Delta
r = s = .1 # da tau2

prob.samp[1,] = gtools::rdirichlet(1, alpha = xi) # pesos
beta.samp[1,,] = matrix(rnorm(p*G,sd=c),ncol=G) # coef reg
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


# -------------------------------------------------------------------------------
#                                    Gibbs Sampling
# -------------------------------------------------------------------------------

phi = 0.3 # "desvio padrao" da proposta de MH
cont = 0 # contador de aceites de MH

library(compiler)
enableJIT(3)

for(i in 2:Q){
 b = -sqrt(nu.samp[i-1]/pi)*(gamma((nu.samp[i-1]-1)/2)/gamma(nu.samp[i-1]/2))
 
 # ----------------------
 # Atualizando alpha
 # ----------------------
 alpha.samp[i] = full_alpha(nu.samp[i-1])
 
 # atualizando pametros por componente
 for(j in 1:G){
  
  # ----------------------
  # Atualizando beta
  # ----------------------
  beta.samp[i,,j] = full_beta(j,tau2.samp[i-1,],Delta.samp[i-1,],
                              u.samp[i-1,],t.samp[i-1,],z.samp[i-1,])
  
  # ----------------------
  # Atualizando tau2
  # ----------------------
  tau2.samp[i,j] = full_tau2(j,beta.samp[i,,],Delta.samp[i-1,],
                             u.samp[i-1,],t.samp[i-1,],z.samp[i-1,])
  
  # ----------------------
  # Atualizando Delta
  # ----------------------
   Delta.samp[i,j] = full_Delta(j,beta.samp[i,,],tau2.samp[i,],
                                u.samp[i-1,],t.samp[i-1,],z.samp[i-1,])
   # aqui tomamos eta = 0 para evitar mais conta, lembrar de 
   # alterar a media na funcao full_delta caso eta != 0.
 }
 
 # ----------------------
 # Atualizando prob
 # ----------------------
 prob.samp[i,] = full_prob(z.samp[i-1,])
 
 # ----------------------
 # Atualizando Z
 # ----------------------
 z.samp[i,] = full_Z(prob.samp[i,],beta.samp[i,,],tau2.samp[i,],Delta.samp[i,],
                     nu.samp[i-1])
  
 # ----------------------
 # Atualizando T
 # ----------------------
 t.samp[i,] = full_T(beta.samp[i,,],tau2.samp[i,],Delta.samp[i,],u.samp[i-1,],
                     z.samp[i,])
 
 # ----------------------
 # Atualizando U
 # ----------------------
 u.samp[i,] = full_U(beta.samp[i,,],tau2.samp[i,],Delta.samp[i,],nu.samp[i-1],
                     t.samp[i,],z.samp[i,])
  
 # ----------------------
 # Atualizando nu
 # ----------------------
 nu.samp[i] = full_nu2(prob.samp[i,],beta.samp[i,,],tau2.samp[i,],Delta.samp[i,],
                       alpha.samp[i],nu.samp[i-1])
 
 cont = ifelse(nu.samp[i]==nu.samp[i-1],cont,cont + 1)
  
 pb$tick(tokens = list(taxa = paste0(formatC(cont/i * 100,2,format="f"),"%")))
};beepr::beep()
