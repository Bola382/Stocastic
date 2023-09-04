# Gera amostras de uma qui-quadrado com nu graus de liberdade
# utilizando soma de normais padrao, requer myrnorm
# n: numero de amostras (natural)
# nu: graus de liberdade (natural)
source("myrnorm.R")

myrchisq = function(n,nu){
 x = matrix(myrnorm(n*nu,0,1)^2, nrow = n, ncol = nu)
 
 rowSums(x)
}
