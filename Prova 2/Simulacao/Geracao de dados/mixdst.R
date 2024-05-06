# densidade misturas de st
# p: vetor 1:G com pesos das comps
# mu, sigma, lambda: vetores 1:G com os respectivos params
# nu: escalar
mixdst = function(x,p,mu,sigma,lambda,nu){
 G = length(p)
 index = 1:G
 
 sum(sapply(index, function(a) p[a] * dst(x, mu[a], sigma[a], lambda[a], nu)))
}

# xx = seq(-10,10, by = .01)
# yy = sapply(xx, function(a) mixdst(a,c(.8,.2),c(-5,5),c(1,1),c(-8,8),3))
# 
# plot(xx,yy,type="l")
