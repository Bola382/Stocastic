# ==================================================================
# Trata o problema de label switching (ta errado)
# as cadeias geradas nao estao disponiveis por conta de seu tamanho
# ===================================================================
nsamp = nrow(prob.samp)

# montando uma matriz de parametros
probvec = c(prob.samp)
betamat = matrix(NA,nrow = nsamp*G, ncol = p) # cada coluna vai ser um dos coefs
                                              # misturando os coefs de comps diferentes
for(j in 1:p){
 betamat[,j] = c(beta.samp[,j,])
}

tau2vec = c(tau2.samp)
Deltavec = c(Delta.samp)

# matriz de parametros
theta = cbind(probvec,betamat,tau2vec,Deltavec)

# novos rotulos
set.seed(2)
aux = kmeans(theta,centers=G,iter.max = 150)
labels = aux$cluster
aux$size # devemos ter (Q-burn)/thin + 1  em cada cluster

# re-rotulando os parametros
for(j in 1:G){
 prob.samp[,j] = probvec[labels==j]
 beta.samp[,,j] = betamat[labels==j,]
 tau2.samp[,j] = tau2vec[labels==j]
 Delta.samp[,j] = Deltavec[labels==j]
}

rm(list=setdiff(ls(), c("prob.samp","beta.samp","tau2.samp","Delta.samp","nu.samp",
                        "alpha.samp","z.samp","u.samp","t.samp","labels","Q","p")))
