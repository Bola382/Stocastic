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

# acumulada st
pst = function(q,mu,sigma,lambda,nu){
 dst2 = function(x){dst(x,mu,sigma,lambda,nu)}
 integrate(dst2, lower = -Inf, upper = q, rel.tol = 1e-10)$value
}

# funcao quantil st
qst = function(p,mu,sigma,lambda,nu){
 aux = function(q){p - pst(q,mu,sigma,lambda,nu)}
 
 uniroot(aux, interval = c(-50,50), extendInt = "yes")$root
 
}

# para fazer o grafico quantil-quantil
qqst = function(sample,mu,sigma,lambda,nu){
 n = length(sample)
 sample = sort(sample)
 teo = sapply((1:n-.5)/n, function(a) qst(a,mu,sigma,lambda,nu))
 
 return(cbind(ord = sample, teo = teo))
}

# ------------------------------------------------------------------
#                      Representacao hierarquica
# ------------------------------------------------------------------

# gera amostras da st 
# usando a representacao hierarquica
source("rstHier.R")

# testando uma amostra de tamanho 50
set.seed(1)
samp_hier = rstHier(50, mu, sigma, lambda, nu) # amostra

qq_hier = qqst(samp_hier,mu,sigma,lambda,nu) # qq

plot(qq_hier[,1],qq_hier[,2]);abline(0,1) # qq

hist(samp_hier, freq=F, breaks=16)
curve(dst(x, mu, sigma, lambda, nu), add = T)

plot(ecdf(samp_hier)) # acumulada empirica
xx = seq(min(samp_hier),max(samp_hier), by = .5) # acumulada teo
pst2 = function(q){pst(q,mu,sigma,lambda,nu)} # acumulada teo
lines(xx, sapply(xx, pst2)) # acumulada teo

ks.test(samp_hier,Vectorize(pst2)) # teste ks

# testando R amostras de tamanho n
R = 100
n = 1000

samps_hier = list() # amostras
times_hier = NULL # tempos de geracao
pvalues_hier = NULL # pvalores de cada amostra
qqs_hier = array(NA, c(R, n, 2), # salva os pontos do qqplot
                 dimnames = list(rep = 1:R, samp = 1:n, c("ord","teo")))

# barra de progresso
pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                     max = R, # Maximum value of the progress bar
                     style = 3,    # Progress bar style (also available style = 1 and style = 2)
                     width = 50,   # Progress bar width. Defaults to getOption("width")
                     char = "=")   # Character used to create the bar

set.seed(4) # seed 3 para n = 100, 4 para n = 1000, 5 para n = 5000
for(i in 1:R){
 t1 = proc.time()[3]
 samps_hier[[i]] = rstHier(n, mu, sigma, lambda, nu)
 times_hier[i] = proc.time()[3] - t1
 
 qqs_hier[i,,] = qqst(samps_hier[[i]], mu, sigma, lambda, nu)
 
 pvalues_hier[i] = ks.test(samps_hier[[i]],Vectorize(pst2))$p.value
 
 setTxtProgressBar(pb, i)
};close(pb);beepr::beep()

summary(times_hier)
plot(pvalues_hier, pch = 16); abline(h=.05, col =2, lty =2)
sum(pvalues_hier < .05)

apply(qqs_hier, 3, min)
apply(qqs_hier, 3, max)

plot(1, type = "n", xlim = c(-100,0), ylim = c(-60, -2),
     xlab = "quantis amostrais", ylab = "quantis teóricos",
     main = "Gráfico quantil-quantil")
apply(qqs_hier, 1, points, pch = 16, cex = .6); abline(0,1,col=2,lty=2)
legend("topleft", legend = "y = x", lty = 2, col = 2)

# integracao de monte carlo
n = 10000

set.seed(10)
samp_MC = rstHier(n, mu, sigma, lambda, nu)

qs = seq(-20,-10,by=.01)

quants = sapply(qs, function(a) mean(samp_MC <= a)); quants

# prob no quantil
quants[which.min(abs(quants-.5))]
qs[which.min(abs(quants-.5))]

qst(.5,mu,sigma,lambda,nu)

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
y = rstRej(100, mu, sigma, lambda, nu, myrnorm,env_rej, -15, 20)
hist(y$samples, freq = F, breaks = 30, add = T)
100/y$interacoes

samp_rej = y$samples

qq_rej = qqst(samp_rej,mu,sigma,lambda,nu)

plot(qq$teoretical,qq$sample);abline(0,1)

hist(samp_hier, freq=F, breaks=16)
curve(dst(x, mu, sigma, lambda, nu), add = T)

ks.test(samp_rej,Vectorize(pst2))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~
#        Com squeeze
# ~~~~~~~~~~~~~~~~~~~~~~~~~~

# vizualizando objetivo, proposta e squeeze
plot(1, xlim = c(-60,10), ylim = c(-.01, .1),type = "n",
     xlab = "x", ylab = "Densidade")
curve(dst(x, mu, sigma, lambda, nu), add = T)

# funcao envelope
env_rej = function(x) 3.88*dnorm(x, -15, 20) 
# funcao squeeze (olhar applet do geogebra)
sqz_rej = function(x){
 if(x > -38 & x < -5.8){
 -0.0000000061715*x^6 - 0.0000009267749*x^5 - 0.0000560775342*x^4-
  0.0017399225511*x^3 - 0.0288821956928*x^2 - 0.2364142400584*x -
  0.6715760552873
 }else{0}
}

curve(env_rej(x), col = 2, lty = 2, add = T)
xx = seq(-60,10,by=.5)
lines(xx,sapply(xx,sqz_rej), col = 2)
legend("topleft", legend = c("Função objetivo", 
                             "Squeeze",
                             "Envelope proposto"), 
       lty = c(1,1,2), col = c(1,2,2))

# verificando funcionamento
source("rstRej_sqz.R")

source("myrnorm.R")
set.seed(1)
y = rstRej_sqz(100, mu, sigma, lambda, nu, myrnorm,env_rej, sqz_rej, -15, 20)
100/y$iteracoes # taxa de aceitacao

samp_rejSqz = y$samples

qq_rej = qqst(samp_rejSqz,mu,sigma,lambda,nu)

plot(qq_rej[,1],qq_rej[,2], 
     xlab = "quantis amostrais", 
     ylab = "quantis teóricos",
     main = "Gráfico quantil-quantil",
     pch=16, cex = .6);abline(0,1)
legend("topleft", legend = "y = x", lty = 1)

plot(ecdf(samp_rejSqz)) # acumulada empirica
xx = seq(min(samp_rejSqz),max(samp_rejSqz), by = .5) # acumulada teo
pst2 = function(q){pst(q,mu,sigma,lambda,nu)} # acumulada teo
lines(xx, sapply(xx, pst2)) # acumulada teo

ks.test(samp_rejSqz,Vectorize(pst2))

# testando R amostras de tamanho n
R = 100
n = 5000

samps_rej = list() # amostras
times_rej = NULL # tempos de geracao
niter_rej = NULL # numero de iteracoes
pvalues_rej = NULL # pvalores de cada amostra
qqs_rej = array(NA, c(R, n, 2), # salva os pontos do qqplot
                 dimnames = list(rep = 1:R, samp = 1:n, c("ord","teo")))

# barra de progresso
pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                     max = R, # Maximum value of the progress bar
                     style = 3,    # Progress bar style (also available style = 1 and style = 2)
                     width = 50,   # Progress bar width. Defaults to getOption("width")
                     char = "=")   # Character used to create the bar

set.seed(4) # trocar seed ao trocar n
for(i in 1:R){
 t1 = proc.time()[3]
 y = rstRej_sqz(n, mu, sigma, lambda, nu, myrnorm,env_rej, sqz_rej, -15, 20)
 
 samps_rej[[i]] = y$samples
 niter_rej[i] = y$iteracoes
 times_rej[i] = proc.time()[3] - t1
 
 qqs_rej[i,,] = qqst(samps_rej[[i]], mu, sigma, lambda, nu)
 
 pvalues_rej[i] = ks.test(samps_rej[[i]],Vectorize(pst2))$p.value
 
 setTxtProgressBar(pb, i)
};close(pb);beepr::beep()

summary(times_rej)
plot(pvalues_rej, pch = 16); abline(h=.05, col =2, lty =2)
sum(pvalues_rej<.05)

summary(n/niter_rej)

apply(qqs_rej, 3, min)
apply(qqs_rej, 3, max)

plot(1, type = "n", xlim = c(-105,1), ylim = c(-70, -1),
     xlab = "quantis amostrais", ylab = "quantis teóricos",
     main = "Gráfico quantil-quantil")
apply(qqs_rej, 1, points, pch = 16, cex = .6); abline(0,1,col=2,lty=2)
legend("topleft", legend = "y = x", lty = 2, col = 2)

