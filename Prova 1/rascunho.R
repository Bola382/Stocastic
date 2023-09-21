curve(dst(x,mu,sigma,lambda,nu), from = -60, to = 10, ylab = "Densidade")
curve(dst(x,mu,sigma,-lambda,nu), from = -10, to = 60, ylab = "Densidade")

h = function(x) log(dst(x,mu,sigma,-lambda,nu))

hlinha = function(x){
 a = (x-mu)/sigma
 aa = a*(-lambda)*sqrt((nu+1)/(a^2+nu))
 
 -((nu+1)/nu*a/sigma)/(1+a^2/nu) + (dt(aa,df=nu+1)*aa*(1/(x-mu) - (a/sigma)/(a^2+nu)))/pt(aa,df=nu+1)
}

curve(h(x), from = -200, to = 10)
curve(hlinha(x), from = -200, to = 10, ylab = "Derivada de h_2")

abline(a = h(-150)-hlinha(-150)*(-150), b = hlinha(-150), col = 2)
abline(a = h(-50)-hlinha(-50)*(-50), b = hlinha(-50), col = 2)
abline(a = h(0)-hlinha(0)*(0), b = hlinha(-0), col = 2)

xvec = c(-150,-100,-50,-10,5)
yvec = h(xvec)

points(xvec,yvec)

un = function(x){
 if(x<=xvec[1]){ # reta entre x1 e x2
  a = (yvec[2]-yvec[1])/(xvec[2]-xvec[1])
  b = yvec[1] - a * xvec[1]
  a*x + b
 }else{
  if(x > xvec[1] & x <= xvec[2]){ # reta entre x2 e x3
   a = (yvec[3]-yvec[2])/(xvec[3]-xvec[2])
   b = yvec[2] - a * xvec[2]
   a*x + b
  }else{
   if(x > xvec[length(xvec)]){# reta entre xn-1 e xn
    a = (yvec[length(xvec)]-yvec[length(xvec)-1])/(xvec[length(xvec)]-xvec[length(xvec)-1])
    b =  yvec[length(xvec)-1] - a * xvec[length(xvec)-1]
    a*x+b
   }else{0}
  }
 }
}
lines(seq(-200,10,by=.5),sapply(seq(-200,10,by=.5),un))

h(-200) - 2*h(-180) + h(-150)
