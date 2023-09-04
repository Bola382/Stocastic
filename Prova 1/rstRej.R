# Gera amostras da st com parametros mu, sigma, lambda, nu 
# usando o algoritmo da rejeicao, requer dst.R
# n: numero de amostras (natural)
# mu, sigma, lambda, nu: parametros (real, real positivo, real, natural)
# dist: funcao que gera amostras da distribuicao proposta
# env: funcao que da o envelope proposto, deve ser proporcional a fdp de dist
# ...: parametros opcionais de dist

rstRej = function(n, mu, sigma, lambda, nu, dist, env, ...){
 inter = 0    # numero de interacoes
 nsamp = 0    # amostras da st geradas
 samp  = NULL # vetor de amostras
 
 while(nsamp < n){
  x = dist(1, ...)
  u = runif(1)
  
  p = dst(x, mu, sigma, lambda, nu)/env(x)
  
  if(p>1) stop("p > 1, verificar envelope.")
  
  if(u <= p){samp = c(samp,x); nsamp = nsamp+1}
  
  inter = inter + 1
 }
 return(list(samples = samp, interacoes = inter))
}
