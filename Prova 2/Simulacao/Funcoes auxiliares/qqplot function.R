# Calcula o residuo quantilico para o modelo de misturas de regressoes t assimetrica
# dados: matriz de dados, coluna 1 respostas e restantes covariaveis
# restante dos argumentos vem das amostras a posteriori
qqplotmix = function(dados,beta,prob,sigma,lambda,nu){
 n = nrow(dados)
 p = ncol(dados)
 Q = length(nu)
 
 y = dados[,1] # respostas
 X = cbind(1,dados[,-1]) # covs com intercepto
 
 G = ncol(beta)/p
 b = sapply(1:Q, function(a) ifelse(nu[a] > 1,-exp(log(nu[a])/2 - log(pi)/2 + lgamma((nu[a]-1)/2) - lgamma(nu[a]/2)),-sqrt(nu[a]/pi)*(gamma((nu[a]-1)/2)/gamma(nu[a]/2))))
 Delta = sigma*lambda/sqrt(1+lambda^2)
 bDelta = matrix(rep(b,G),Q,G)*Delta
 
 qq = matrix(NA, nrow = Q, ncol = n)
 
 for(iter in 1:Q){
  # parametro de localizacao para cada i na iteracao iter
  mean_vec = X %*% matrix(as.numeric(beta[iter,]),p,G) + matrix(as.numeric(bDelta[iter,]), n, G, byrow = T)# n x G
  
  # residuo quantilico
  qq[iter,] = sort(qnorm(sapply(1:n, function(a) mixpst(y[a], as.numeric(prob[iter,]), mean_vec[a,], as.numeric(sigma[iter,]), as.numeric(lambda[iter,]), nu[iter]))))
 }
 tab = cbind(colMeans(qq),t(apply(qq, 2, function(a) quantile(a, c(.025,.975)))))
 colnames(tab) = c("res", "LCI","UCI")
 
 return(tab)
}