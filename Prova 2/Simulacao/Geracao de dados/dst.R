# densidade da st por azallini
# mesmos parametros de rstHier

dst = function(x,mu,sigma,lambda,nu){
 a = (x-mu)/sigma
 2*dt(a, nu)/sigma*pt(a*lambda*sqrt((nu+1)/(a^2+nu)), nu+1)
}