# gera amostras de uma t assimetrica utilizando a representacao hierarquica

# n: numero de amostras (natural)
# mu: parametro de locacao (real)
# sigma: raiz do parametro de escala (real positivo)
# lambda: parametro de forma (real)
# nu: graus de liberdade (natural por limitacoes)

rstHier = function(n, mu, sigma, lambda, nu){
 delta = lambda/sqrt(1+lambda^2)
 Delta = sigma*delta
 tau = sigma*sqrt(1-delta^2)
 
 # calcula a transformacao de Box-Muller
 # a: vetor bidimensional
 boxMuller = function(a){
  c(sqrt(-2*log(a[1]))*cos(2*pi*a[2]), sqrt(-2*log(a[1]))*sin(2*pi*a[2]))
 }
 
 # ------------------------------------
 # Gerando valores de U
 # ------------------------------------
 
 # gerando de uniformes para gerar amostras
 # de U,
 # se n*nu for impar n1 = n*nu+1
 # cc n1 = n*nu
 # assim sempre temos pares
 n1 = ifelse((n*nu)%%2,n*nu+1,n*nu)
 
 unifs = matrix(runif(n1), nrow = n1/2, 2)
 
 # aplica por linhas de unifs a tranformacao de Box-Muller
 norms1 = c(apply(unifs, 1, boxMuller))
 
 # soma de normais ao quadrado
 # u ~ gama(nu/2,nu/2)
 u = rowSums(matrix(norms1[1:(n*nu)]^2, nrow=n))/nu
 
 # ------------------------------------
 # Gerando valores de T0 e T1
 # ------------------------------------
 
 unifs2 = matrix(runif(2*n), nrow = n, 2)
 
 # aplica por linhas de unifs2 a tranformacao de Box-Muller
 norms2 = t(apply(unifs2, 1, boxMuller))
 
 t0 = abs(norms2[,1])
 t1 = norms2[,2]
  
 mu + (Delta*t0 + tau*t1)/sqrt(u)
}