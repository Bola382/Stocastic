# gera amostras de uma t assimetrica utilizando a representacao hierarquica
# requer funcao myrnorm e myrchisq
# n: numero de amostras (natural)
# mu: parametro de locacao (real)
# sigma: raiz do parametro de escala (real positivo)
# lambda: parametro de forma (real)
# nu: graus de liberdade (natural por limitacoes do myrchisq)
source("myrnorm.R")
source("myrchisq.R")

rstHier = function(n, mu, sigma, lambda, nu){
 delta = lambda/sqrt(1+lambda^2)
 Delta = sigma*delta
 tau = sigma*sqrt(1-delta^2)
 
 u = myrchisq(n, nu)/nu # gama(nu/2,nu/2)
 x = myrnorm(n, 0, 1) # normal(0,1)
 t = abs(myrnorm(n, 0, 1)) # normal(0,1) truncada em 0
 
 mu + (Delta*t + tau*x)/sqrt(u)
}
