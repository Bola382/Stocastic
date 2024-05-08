# ==============================================================================
# Gerando amostras da posteriori utilizando GS
# salva resultados em um .txt
# ==============================================================================
# necessita carregar a pasta "Funcoes auxiliares" e a funcao "dst.R"
# e do pacote "compiler"

# folder: nome da pasta em qual serao salvos os resultados de cada replicacao
# repl: numero da replicacao atual
# data: matriz de dados, primeira coluna respostas, demais valores das covs
# G: numero de componentes
# phi: "desvio padrao" da proposta de MH para nu
# Q: numero de iteracoes do algoritmo
# burn: amostras a serem descartadas
# thin: tamanho dos saltos pos burnin

gibbs = function(folder,repl,dados,G,phi=.3,Q,burn,thin){
 n = nrow(dados)
 p = ncol(dados) # numero de betas
 
 y = dados[,1] # resp
 X = as.matrix(cbind(1,dados[,-1])) # covariaveis com intercepto
 
 # ~~~~~~~~~~~~~~~~
 # valores iniciais
 # ~~~~~~~~~~~~~~~~
 
 # para guardar amostras
 prob.samp = tau2.samp = Delta.samp = eta.samp = matrix(NA, nrow = 2, ncol = G)
 u.samp = t.samp = z.samp = matrix(NA, nrow = 2, ncol = n)
 beta.samp = array(NA, dim = c(2,p,G), dimnames = list(c("old","new"), 1:p, 1:G))
 
 # hiperparametros
 xi = rep(1,G) # prob
 c = 10 # da beta (c <- sqrt(c) das minhas contas)
 omega = 10 # da Delta
 sd.eta = 10 # dp do hiperparametro da Delta
 r = s = .1 # da tau2
 
 prob.samp[1,1:G] = gtools::rdirichlet(1, alpha = xi) # pesos
 beta.samp[1,,1:G] = matrix(rnorm(p*G,sd=c),ncol=G) # coef reg
 tau2.samp[1,1:G] = 1/rgamma(G,r,s)# escala
 eta.samp[1,1:G] = rnorm(G, 0, sd.eta) # hiperparametro da priori de Delta
 Delta.samp[1,1:G] = rnorm(G,eta,omega)# forma
 alpha.samp = runif(1,.02,.5)# hiperparametro da priori de nu
 nu.samp = rexp(1,rate=alpha.samp)# gl
 
 
 etaFull = sd.eta^2/(sd.eta^2+omega^2) # termo que aparece na media da condicional completa de eta
 etaVar = etaFull*omega^2 # variancia da cond completa de eta
  
 # centralizacao da media
 b = -sqrt(nu.samp/pi)*(gamma((nu.samp-1)/2)/gamma(nu.samp/2))
 
 u.samp[1,] = rgamma(n, shape = nu.samp/2, rate = nu.samp/2) # t e u da representacao aumentada
 t.samp[1,] = abs(rnorm(n))/sqrt(u.samp[1,])
 z.samp[1,] = sample(1:G,n,prob=prob.samp[1,1:G], replace = T) # latente que indica os grupos
 
 # ~~~~~~~~~~~~~~~~~~~
 # Config output
 # ~~~~~~~~~~~~~~~~~~~
 cat(paste0(rep(paste0("beta",1:p),G),"_",rep(1:G,each=p)),"\n", file = paste0("Outputs/Paralelo/",file,"/",repl,"/beta.txt"))
 cat(paste0("tau2_",1:G),"\n", file = paste0("Outputs/Paralelo/",file,"/",repl,"/tau2.txt"))
 cat(paste0("Delta_",1:G),"\n", file = paste0("Outputs/Paralelo/",file,"/",repl,"/Delta.txt"))
 cat(paste0("p_",1:G),"\n", file = paste0("Outputs/Paralelo/",file,"/",repl,"/prob.txt"))
 cat("nu","alpha",paste0("eta_",1:G),"\n", file = paste0("Outputs/Paralelo/",file,"/",repl,"/param.txt"))
 cat(paste0("z_",1:n),"\n", file = paste0("Outputs/Paralelo/",file,"/",repl,"/z.txt"))
 cat("u.bar","t.bar","\n", file = paste0("Outputs/Paralelo/",file,"/",repl,"/ut.txt"))
 cat("time","rate","\n", file = paste0("Outputs/Paralelo/",file,"/",repl,"/monitor.txt"))
 # -------------------------------------------------------------------------------
 #                                    Amostrador
 # -------------------------------------------------------------------------------
 
 enableJIT(3)
 t_tmp = Sys.time()
 cont = 0
 for(i in 2:Q){
  U = diag(u.samp[1,])
  m = NULL
  #m = table(z.samp[i-1,]) # n de amostras por comp
  
  # ----------------------
  # Atualizando alpha
  # ----------------------
  
  # gama truncada em (0.02,0.5)
  alpha.samp[2] = Runuran::urgamma(1, shape = 2, scale = 1/nu.samp[1], lb = .02, ub = .5)
  
  #print(u.samp[i-1,1:10])
  
  # ----------------------
  # Atualizando eta
  # ----------------------
  
  eta.samp[2,] = rnorm(G,etaFull*Delta.samp[1,], etaVar)
  
  # atualizando pametros por componente
  for(j in 1:G){
   index = which(z.samp[1,]==j)
   mj = length(index)
   m[j] = mj
   Xj = X[index,] # dim = mjXp
   Uj = U[index,index] # dim = mjXmj
   yj = y[index] # dim = mjX1
   tj = t.samp[1, index] # dim = mjX1
   
   if(mj==1){
    XU = as.matrix(Xj)%*%Uj # t(Xj)%*%Uj
   }else{
    XU = crossprod(Xj,Uj) # t(Xj)%*%Uj
   } 
   aux_y = yj-(b+tj)*Delta.samp[1,j] 
   sigma_inv = XU%*%Xj + diag(tau2.samp[1,j]/c^2,nrow = p)
   
   # ----------------------
   # Atualizando beta
   # ----------------------
   s_sigma = chol2inv(chol(sigma_inv)) #solve(sigma_inv)
   mmm = solve(sigma_inv,XU%*%aux_y) 
   ppp = tau2.samp[1,j]*s_sigma#solve(sigma_inv)
   # [iter,coefs,componente]
   beta.samp[2,,j] = MASS::mvrnorm(1, mu = mmm, 
                                   Sigma = ppp)
   #if(j==1){cat(i,"\t", beta.samp[i,,j], "\t med",mmm[1],"\t var", ppp[1,1], "\n")}
   # ----------------------
   # Atualizando tau2
   # ----------------------
   S3 = sum(u.samp[1,index]*(yj-Xj%*%beta.samp[2,,j]-(b+tj)*Delta.samp[1,j])^2)
   
   tau2.samp[2,j] = LearnBayes::rigamma(1, mj/2+r, (S3+s)/2)
   #MCMCpack::rinvgamma(1, mj/2+r, 2/(S3+s))
   #1/rgamma(1,shape = m[j]/2+r, rate = (S3 + s)/2)
   
   # ----------------------
   # Atualizando Delta
   # ----------------------
   S1 = sum(u.samp[1,index]*(yj-Xj%*%beta.samp[2,,j])*(b+tj))
   S2 = sum(u.samp[1,index]*(b+tj)^2)
   denom = omega^2*S2 + tau2.samp[2,j]
   
   Delta.samp[2,j] = rnorm(1, (omega^2*S1 + eta.samp[2,j]*tau2.samp[2,j])/denom, omega*sqrt(tau2.samp[2,j]/denom))
  }
  
  # ----------------------
  # Atualizando prob
  # ----------------------
  
  prob.samp[2,] = gtools::rdirichlet(1, alpha = xi + m)
  
  # ------------------------------
  # atualizando variaveis latentes
  # ------------------------------
  aux_mu = X%*%beta.samp[2,,]+matrix(rep(b*Delta.samp[2,],n), nrow = n, ncol = G, byrow = T)
  aux_sig = sqrt(tau2.samp[2,]+Delta.samp[2,]^2)
  aux_lam = Delta.samp[2,]/sqrt(tau2.samp[2,])
  
  for(k in 1:n){
   # ----------------------
   # Atualizando Z
   # ----------------------
   aux_prob = prob.samp[2,]*sapply(1:G,function(a) dst(y[k], aux_mu[k,a], aux_sig[a], aux_lam[a], nu.samp[1]))
   z.samp[2,k] = sample(1:G,1, prob = aux_prob) # prob eh normalizado internamente
   
   # ----------------------
   # Atualizando T
   # ----------------------
   aux_yxbeta = y[k]-aux_mu[k,z.samp[2,k]]
   aux_denom = Delta.samp[2,z.samp[2,k]]^2 + tau2.samp[2,z.samp[2,k]]
   aux_mean = (aux_yxbeta*Delta.samp[2,z.samp[2,k]])/aux_denom
   aux_sd = sqrt(tau2.samp[2,z.samp[2,k]]/(u.samp[1,k]*aux_denom))
   # media antes do truncamento muito distante do lb
   
   aux_truncnorm = rtruncnorm(1,-aux_mean/aux_sd)$amostra
   
   t.samp[2,k] = aux_mean+aux_sd*aux_truncnorm
   # ----------------------
   # Atualizando U
   # ----------------------
   D1 =(aux_yxbeta-t.samp[2,k]*Delta.samp[2,z.samp[2,k]])^2/(2*tau2.samp[2,z.samp[2,k]]) 
   u.samp[2,k] = rgamma(1, shape = nu.samp[1]/2+1, rate = D1 + (t.samp[2,k]^2+nu.samp[1])/2)
  }
  # ----------------------
  # Atualizando nu
  # ----------------------
  
  # Passo de MH
  prop = rlnorm(1, log(nu.samp[1]), phi)
  
  aceit = min(1,exp(fullnu(prop,alpha.samp[2],prob.samp[2,],y,aux_mu,aux_sig,aux_lam) + log(prop) - fullnu(nu.samp[1],alpha.samp[2],prob.samp[2,],y,aux_mu,aux_sig,aux_lam) - log(nu.samp[1])))
  if(aceit > 1 | aceit <0){stop("Problema no passo MH")}
  if(runif(1)<=aceit){
   nu.samp[2] = prop
   cont = cont + 1
  }else{
   nu.samp[2] = nu.samp[1]
  }
  
  # centralizacao da media
  b = ifelse(nu.samp[2] > 1,-exp(log(nu.samp[2])/2 - log(pi)/2 + lgamma((nu.samp[2]-1)/2) - lgamma(nu.samp[2]/2)),-sqrt(nu.samp[2]/pi)*(gamma((nu.samp[2]-1)/2)/gamma(nu.samp[2]/2)))
  
  # output
  if(i>burn & i%%thin==0){
   cat(c(beta.samp[2,,]),"\n", file = paste0("Outputs/Paralelo/",file,"/",repl,"/beta.txt"), append=T)
   cat(tau2.samp[2,],"\n", file = paste0("Outputs/Paralelo/",file,"/",repl,"/tau2.txt"), append=T)
   cat(Delta.samp[2,],"\n", file = paste0("Outputs/Paralelo/",file,"/",repl,"/Delta.txt"), append=T)
   cat(prob.samp[2,],"\n", file = paste0("Outputs/Paralelo/",file,"/",repl,"/prob.txt"), append=T)
   cat(nu.samp[2],alpha.samp[2],eta.samp[2,],"\n", file = paste0("Outputs/Paralelo/",file,"/",repl,"/param.txt"), append=T)
   cat(z.samp[2,],"\n", file = paste0("Outputs/Paralelo/",file,"/",repl,"/z.txt"), append=T)
   cat(colMeans(cbind(u.samp[2,],t.samp[2,])),"\n", file = paste0("Outputs/Paralelo/",file,"/",repl,"/ut.txt"), append=T)
  }
  # o novo se torna antigo na proxima iteracao
  z.samp[1,] = z.samp[2,]
  beta.samp[1,,] = beta.samp[2,,]
  tau2.samp[1,] = tau2.samp[2,]
  Delta.samp[1,] = Delta.samp[2,]
  u.samp[1,] = u.samp[2,]
  t.samp[1,] = t.samp[2,]
  alpha.samp[1] = alpha.samp[2]
  prob.samp[1,] = prob.samp[2,]
  nu.samp[1] = nu.samp[2]
  
  # limpando o novo
  z.samp[2,] = u.samp[2,] = t.samp[2,] = rep(NA,n)
  prob.samp[2,] = tau2.samp[2,] = Delta.samp[2,] = rep(NA, G)
  alpha.samp[2] = nu.samp[2] = NA
  beta.samp[2,,] = matrix(NA,p,G)
  
 };time_ok = Sys.time()-t_tmp
 cat(c(time_ok,cont/Q),"\n", file = paste0("Outputs/Paralelo/",file,"/",repl,"/monitor.txt"), append=T)
}
