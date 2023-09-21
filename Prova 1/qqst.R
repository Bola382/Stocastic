qqst = function(sample,mu,sigma,lambda,nu){
 n = length(sample)
 sample = sort(sample)
 teo = sapply(1:n/n, function(a) qst(a,mu,sigma,lambda,nu))
 
 return(sample = sample, teoretical = teo)
}