# Gera amostras da st com parametros mu, sigma, lambda, nu 
# usando o algoritmo da rejeicao, requer dst.R
# n: numero de amostras (natural)
# mu, sigma, lambda, nu: parametros (real, real positivo, real, natural)
# dist: funcao que gera amostras da distribuicao proposta
# o primeiro argumento de dist deve ser o numero de amostras
# env: funcao que da o envelope proposto, deve ser proporcional a fdp de dist
# sqz: funcao squeeze
# ...: parametros opcionais de dist

rstRej_sqz = function(n, mu, sigma, lambda, nu, dist, env, sqz, ...){
 inter = 0    # numero de interacoes
 nsamp = 0    # amostras da st geradas
 samp  = NULL # vetor de amostras
 
 while(nsamp < n){
  x = dist(1, ...)
  u = runif(1)
  
  p = sqz(x)/env(x)
  
  if(p>1 | p < 0) stop("p nao esta no (0,1), verificar envelope e squeeze.")
  
  if(u <= p){
   samp = c(samp,x); nsamp = nsamp+1
  }else{
   p = dst(x, mu, sigma, lambda, nu)/env(x)
   
   if(p>1 | p < 0) stop("p nao esta no (0,1), verificar envelope.")
   if(u <= p){samp = c(samp,x); nsamp = nsamp+1}
  }
  
  inter = inter + 1
 }
 return(list(samples = samp, iteracoes = inter))
}
