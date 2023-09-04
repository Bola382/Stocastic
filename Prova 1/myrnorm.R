# gera amostras de uma distribuicao normal utilizando Box e Muller
# n: numero de amostras (n natural)
# mu: media da normal (mu real)
# sigma: desvio padrao da normal (sigma real positivo)

myrnorm = function(n, mu, sigma){
 # se n for impar n1 = n+1
 # cc n1 = n
 # assim sempre temos pares de u_i, u_j
 n1 = ifelse(n%%2,n+1,n)
 u = matrix(runif(n1), nrow = n1/2, 2)
 
 y = mu + sigma*c(apply(u, 1, function(a) c(sqrt(-2*log(a[1]))*cos(2*pi*a[2]), sqrt(-2*log(a[1]))*sin(2*pi*a[2]))))
 
 return(y[1:n])
}
