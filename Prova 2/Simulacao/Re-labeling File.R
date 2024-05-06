# ==================================================================
# Trata o problema de label switching (Sylvia Fruhwirth-Schnatter)
# ===================================================================
# file: nome da pasta em Outputs/Paralelo onde estao salvos os resultados
# repl: qual replicacao escolher
# try: numero de tentativas maximas do kmeans

relabel = function(file,repl,try=10){
  caminho = paste0("Outputs/Paralelo/",file,"/",repl,"/")
  caminho2 = paste0("Outputs/Paralelo/",file,"/",repl,"/corrigido/")
  dir.create(caminho2,showWarnings = F)
  
  beta.samp = read.table(paste0(caminho,"beta.txt"), h =T)
  tau2.samp = read.table(paste0(caminho,"tau2.txt"), h =T)
  Delta.samp = read.table(paste0(caminho,"Delta.txt"), h =T)
  prob.samp = read.table(paste0(caminho,"prob.txt"), h =T)
  z.samp = read.table(paste0(caminho,"z.txt"), h =T)
  ut = read.table(paste0(caminho,"ut.txt"), h =T)
  umean.samp = ut$u.bar; tmean.samp = ut$t.bar
  param = read.table(paste0(caminho,"param.txt"), h =T)
  nu.samp = param$nu; alpha.samp = param$alpha
  
  # organizando beta.samp em um array
  Q = nrow(beta.samp)
  G = ncol(tau2.samp)
  n = ncol(z.samp)
  p = ncol(beta.samp)/G # numero de betas (covs + 1)
  beta.aux = array(NA, dim = c(Q,p,G), dimnames = list(1:Q,1:p,1:G))
  beta.samp = as.matrix(beta.samp)
  enableJIT(3)
  for(i in 1:Q){
    qual = c(1:(p-1),0)
    for(j in 1:p){
      quais = which(1:ncol(beta.samp) %% p == qual[j])
      beta.aux[i,j,] = beta.samp[i,quais]
    }
  }
  beta.samp = beta.aux; rm("beta.aux") 
  
  # ==================================================================
  #                       Tratando label-switching
  # ==================================================================
  
  # montando uma matriz de parametros
  nsamp = length(nu.samp)
  
  betamat = matrix(NA,nrow = nsamp*G, ncol = p) # cada coluna vai ser um dos coefs
  # misturando os coefs de comps diferentes, ou seja, coluna 1 beta1, qual comp? sim.
  for(j in 1:p){
    betamat[,j] = c(beta.samp[,j,1:G])
  }
  
  tau2vec = c(as.matrix(tau2.samp[,1:G]))
  Deltavec = c(as.matrix(Delta.samp[,1:G]))
  # ficam organizados por iteracao X cluster
  
  # matriz de parametros
  theta = scale(cbind(betamat,tau2vec,Deltavec))
  
  # novos rotulos
  aux = kmeans(theta,centers=G,iter.max = 5000)
  labels = aux$cluster
  new_label = matrix(labels,nrow=nsamp,ncol=G)
  
  # verificando quais rotulos sao uma permutacao valida de {1,...,G}
  idpermu = unlist(sapply(1:nsamp, function(a) if(identical(sort(new_label[a,]),1:G)){a}))
  npermu = length(idpermu)
  taxa_permu = npermu/nsamp # das amostras com G comps, quantas tem permutacoes validas
  if(taxa_permu <= .5){
    cc = 0
    repeat{
      aux = kmeans(theta,centers=G,iter.max = 5000)
      labels = aux$cluster
      new_label = matrix(labels,nrow=nsamp,ncol=G)
      
      # verificando quais rotulos sao uma permutacao valida de {1,...,G}
      idpermu = unlist(sapply(1:nsamp, function(a) if(identical(sort(new_label[a,]),1:G)){a}))
      npermu = length(idpermu)
      taxa_permu = npermu/nsamp
      cc = cc+1
      if(taxa_permu > .8){break}
      if(taxa_permu <= .8 & cc == try & npermu != 0){break}
      if(taxa_permu <= .8 & cc == try & npermu == 0){
       # salvando
       write.table(colnames(beta.samp), file = paste0(caminho2,"beta.txt"),row.names = F)
       write.table(colnames(Delta.samp), file = paste0(caminho2,"Delta.txt"),row.names = F)
       write.table(colnames(param), file = paste0(caminho2,"param.txt"),row.names = F)
       write.table(colnames(prob.samp), file = paste0(caminho2,"prob.txt"),row.names = F)
       write.table(colnames(tau2.samp), file = paste0(caminho2,"tau2.txt"),row.names = F)
       write.table(colnames(ut), file = paste0(caminho2,"ut.txt"),row.names = F)
       write.table(colnames(z.samp), file = paste0(caminho2,"z.txt"),row.names = F)
       
        return(taxa_permu)
       }
    }
  }
   
  # ajustando os rotulos dos parametros em cada permutacao valida
  
  prob_ok = matrix(NA, nrow = npermu, ncol = G, dimnames = list(1:npermu,paste0("p_",1:G)))
  beta_ok = array(NA, dim = c(npermu,p,G), 
                  dimnames = list(1:npermu,1:p,1:G))
  tau2_ok = matrix(NA, nrow = npermu, ncol = G, dimnames = list(1:npermu,paste0("comp",1:G)))
  Delta_ok = matrix(NA, nrow = npermu, ncol = G, dimnames = list(1:npermu,paste0("comp",1:G)))
  z_ok = matrix(NA, nrow = npermu, ncol = n, dimnames = list(1:npermu,1:n))
  t_ok = tmean.samp[idpermu]
  u_ok = umean.samp[idpermu]
  nu_ok = nu.samp[idpermu]
  
  alpha_ok = alpha.samp[idpermu]
  param = cbind(nu_ok, alpha_ok); colnames(param) = c("nu", "alpha")
  ut_ok = cbind(u_ok,t_ok); colnames(ut_ok) = c("u.bar","t.bar")
  
  for(i in 1:npermu){
    for(j in 1:G){
      prob_ok[i,new_label[idpermu[i],j]] = prob.samp[idpermu[i],j]
      tau2_ok[i,new_label[idpermu[i],j]] = tau2.samp[idpermu[i],j]
      Delta_ok[i,new_label[idpermu[i],j]] = Delta.samp[idpermu[i],j]
      for(k in 1:p){
        beta_ok[i,k,new_label[idpermu[i],j]] = beta.samp[idpermu[i],k,j]
      }
      for(l in 1:n){
        z_ok[i,l] = new_label[idpermu[i],z.samp[idpermu[i],l]]
      }
    }
  }
  
  # convertendo beta_ok em uma matriz npermu x (p*G)
  beta.aux = matrix(NA, nrow = npermu, ncol = p*G, dimnames = list(1:npermu,paste0(rep(paste0("beta",1:p),G),"_",rep(1:G,each=p))))
  
  for(ii in 1:npermu){
    beta.aux[ii,] = c(beta_ok[ii,,])
  }
  
  # salvando
  write.table(beta.aux, file = paste0(caminho2,"beta.txt"),row.names = F)
  write.table(Delta_ok, file = paste0(caminho2,"Delta.txt"),row.names = F)
  write.table(param, file = paste0(caminho2,"param.txt"),row.names = F)
  write.table(prob_ok, file = paste0(caminho2,"prob.txt"),row.names = F)
  write.table(tau2_ok, file = paste0(caminho2,"tau2.txt"),row.names = F)
  write.table(ut_ok, file = paste0(caminho2,"ut.txt"),row.names = F)
  write.table(z_ok, file = paste0(caminho2,"z.txt"),row.names = F)

  return(taxa_permu)
}
