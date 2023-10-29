# ==============================================================================
# Gerando amostras da posteriori utilizando Gibbs
# ==============================================================================

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list=ls())

sapply(list.files("Funcoes auxiliares",pattern="*.R$",full.names=TRUE, 
                  ignore.case=TRUE),source,.GlobalEnv)
source("Geracao de dados/dst.R")
load("Geracao de dados/dados.Rdata")
rm("Phinu")

head(data)

n = nrow(data)
p = ncol(data)-1 # -1 pois a primeira coluna contem os grupos
#obs: p = k+1

y = data[,2] # resp
X = as.matrix(cbind(1,data[,-(1:2)])) # covariaveis com intercepto

# ------------------------
# Inicio do algoritmo
# ------------------------

Q = 50000 # numero de iteracoes de Gibbs
cont = 0 # contador de aceites de MH
phi = 0.5 # parametro de taxa da proposta de MH

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

# ~~~~~~~~~~~~~~~~
# Valores iniciais
# ~~~~~~~~~~~~~~~~

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

# ~~~~~~~~~~~~~~~~~~~~~
# Nomes dos parametros
# ~~~~~~~~~~~~~~~~~~~~~
tmp_paramnames = c("prob","beta","tau2","Delta","alpha","nu","u","t","z")

# ~~~~~~~~~~~~~~~~~~~
# Barra de progresso
# ~~~~~~~~~~~~~~~~~~~
library(progress)
format = "(:spin) [:bar] :percent [Decorrido: :elapsedfull || Estimado: :eta] Taxa de aceitacao :taxa"
pb = progress_bar$new(format, clear = FALSE, total = Q, complete = "=", incomplete = "-", width = 100)

library(compiler)
enableJIT(3)

for(i in 2:Q){
 b = -sqrt(nu.samp[i-1]/pi)*(gamma((nu.samp[i-1]-1)/2)/gamma(nu.samp[i-1]/2))
 
 # valores na ultima iteracao
 tmp_param = list(prob = prob.samp[i-1,],
                  beta = beta.samp[i-1,,],
                  tau2 = tau2.samp[i-1,],
                  Delta = Delta.samp[i-1,],
                  alpha = alpha.samp[i-1],
                  nu = nu.samp[i-1],
                  u = u.samp[i-1,],
                  t = t.samp[i-1,],
                  z = z.samp[i-1,])
 
 # --------------------------------------------------
 # selecionando aleatoriamente a ordem de atualizacao
 # --------------------------------------------------
 
 # 9 e o numero de tipos de parametros, permanece fixo nao importa G
 tmp_ord = replicate(9,sample(1:G))
 tmp_arg = sample(tmp_paramnames)
 colnames(tmp_ord) = tmp_arg 
 
 # para cada tipo de parametro
 for(j in 1:9){
  # atualiza os parametros que nao dependem da componente
  tmp_param = updateGibbs1(tmp_param,tmp_arg[j])
  b = ifelse(nu.samp[i-1]==tmp_param$nu,b,-sqrt(tmp_param$nu/pi)*(gamma((tmp_param$nu-1)/2)/gamma(tmp_param$nu/2)))
  # para cada componente
  for(k in 1:G){
   # atualiza os parametros que dependem da componente
   if(all(tmp_arg[j] != tmp_paramnames[2:4])){break}
   # alteramos o calculo da condicional completa de Delta assumindo que
   # eta = 0 sempre, caso seja necessario alterar lembrar de alterar no arquivo
   # "condicionais completas.R"
   tmp_param = updateGibbs2(tmp_param,tmp_arg[j],tmp_ord[k,tmp_arg[j]])
  }
 }
 
 # salvando valores atualizados
 prob.samp[i,] = tmp_param$prob
 beta.samp[i,,] = tmp_param$beta
 tau2.samp[i,] = tmp_param$tau2
 Delta.samp[i,] = tmp_param$Delta
 alpha.samp[i] = tmp_param$alpha
 nu.samp[i] = tmp_param$nu
 u.samp[i,] = tmp_param$u
 t.samp[i,] = tmp_param$t
 z.samp[i,] = tmp_param$z
 
 cont = ifelse(nu.samp[i]==nu.samp[i-1],cont,cont+1)
 
 pb$tick(tokens = list(taxa = paste0(formatC(cont/i * 100,2,format="f"),"%")))
};beepr::beep()

# library(coda)
# 
# # [iter,coefs,componente]
# traceplot(mcmc(beta.samp[,1,1]))
# acf((beta.samp[seq(3000,Q,100),1,1]))
# mean(beta.samp[seq(3000,Q,100),1,1])
# 
# traceplot(mcmc(nu.samp[9000:Q]))
# acf(nu.samp)
# 
# traceplot(mcmc(Delta.samp[,2]))
# acf((Delta.samp[200:Q,1]), lag.max = 100)