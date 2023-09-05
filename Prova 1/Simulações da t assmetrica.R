# ==================================================================
#           Simulacao de amostras de uma t assimetrica
# ==================================================================
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list=ls())

# selecionando parametros para a geracao
set.seed(1)
mu = rnorm(1, 0, 10); mu
sigma = rgamma(1,10,1); sigma
lambda = rnorm(1, 0, 10); lambda
nu = sample(1:20, 1); nu

# densidade st
source("dst.R")

# ------------------------------------------------------------------
#                      Representacao hierarquica
# ------------------------------------------------------------------

# gera amostras da st 
# usando a representacao hierarquica
source("rstHier.R")

# verificando
samp_hier = rstHier(10000, mu, sigma, lambda, nu)

# ------------------------------------------------------------------
#                        Algoritmo de rejeicao
# ------------------------------------------------------------------

# ~~~~~~~~~~~~~~~~~~~~~~~~~~
#        Sem squeeze
# ~~~~~~~~~~~~~~~~~~~~~~~~~~

# vizualizando objetivo e proposta
plot(1, xlim = c(-60,10), ylim = c(0, .1),type = "n",
     xlab = "x", ylab = "Densidade")
curve(dst(x, mu, sigma, lambda, nu), add = T)
curve(dnorm(x,-15, 20), col = 2, add = T)

env_rej = function(x) 3.88*dnorm(x, -15, 20) # funcao envelope

curve(env_rej(x), col = 2, lty = 2, add = T)
legend("topleft", 
       legend = c("Função objetivo", "Densidade proposta","Envelope proposto"), 
       lty = c(1,1,2), col = c(1,2,2))

# verificando funcionamento
source("rstRej.R")

set.seed(1)
y = rstRej(10000, mu, sigma, lambda, nu, myrnorm,env_rej, -15, 20)
hist(y$samples, freq = F, breaks = 30, add = T)
10000/y$interacoes

# ~~~~~~~~~~~~~~~~~~~~~~~~~~
#        Com squeeze
# ~~~~~~~~~~~~~~~~~~~~~~~~~~

# vizualizando objetivo, proposta e squeeze
plot(1, xlim = c(-60,10), ylim = c(0, .1),type = "n",
     xlab = "x", ylab = "Densidade")
curve(dst(x, mu, sigma, lambda, nu), add = T)
curve(dnorm(x,-15, 20), col = 2, add = T)

# funcao envelope
env_rej = function(x) 3.88*dnorm(x, -15, 20) 
# funcao squeeze (olhar applet do geogebra)
sqz_rej = function(x){
 -0.0000000061715*x^6 - 0.0000009267749*x^5 - 0.0000560775342*x^4-
  0.0017399225511*x^3 - 0.0288821956928*x^2 - 0.2364142400584*x -
  0.6715760552873
}

curve(env_rej(x), col = 2, lty = 2, add = T)
curve(sqz_rej(x), col = 3, lty = 2, add = T)
legend("topleft", legend = c("Função objetivo", 
                             "Densidade proposta",
                             "Envelope proposto"), 
       lty = c(1,1,2), col = c(1,2,2))

# verificando funcionamento
source("rstRej_sqz.R")

set.seed(1)
y = rstRej_sqz(10000, mu, sigma, lambda, nu, myrnorm,env_rej, sqz_rej, -15, 20)
hist(y$samples, freq = F, breaks = 30, add = T)
10000/y$interacoes
