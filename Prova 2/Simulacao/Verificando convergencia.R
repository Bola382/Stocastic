# ==================================================================
# Verifica a convergência de uma cadeia MCMC
# as cadeias geradas nao estao disponiveis por conta de seu tamanho
# ==================================================================
library(coda)

rm(list=setdiff(ls(), c("prob.samp","beta.samp","tau2.samp","Delta.samp","nu.samp",
                        "alpha.samp","z.samp","u.samp","t.samp","G","Q","p")))

burn = Q/2 # numero de amostras descartadas
thin = 10 # consideramos de thin em thin amostras

index = seq(burn,Q,thin)

prob.samp = prob.samp[index,]
beta.samp = beta.samp[index,,]
tau2.samp = tau2.samp[index,]
Delta.samp = Delta.samp[index,]
nu.samp = nu.samp[index]
alpha.samp = alpha.samp[index]

z.samp = z.samp[index,]
t.samp = t.samp[index,]
u.samp = u.samp[index,]

# --------------------------------------------
# prob
# --------------------------------------------

par(mfrow = c(1,G), mar = c(5.1, 4.1, 4.1, 2.1))
for(j in 1:G){
 traceplot(mcmc(prob.samp[,j]), xlab = "Iterações",
           main = eval(expression(paste("p",j))))
}
for(j in 1:G){
 acf(prob.samp[,j], xlab = "Iterações",
           main = eval(expression(paste("p",j))))
}

# --------------------------------------------
# Beta
# --------------------------------------------

par(mfrow = c(p,G), mar = c(4,4,3,1))
for(i in 1:p){
 for(j in 1:G){
  traceplot(mcmc(beta.samp[,i,j]), xlab="Iterações", 
            main = eval(expression(paste("beta",i,j))))
 }
}

# --------------------------------------------
# tau2
# --------------------------------------------

par(mfrow = c(1,G), mar = c(5.1, 4.1, 4.1, 2.1))
for(j in 1:G){
 traceplot(mcmc(tau2.samp[1000:Q,j]), xlab = "Iterações",
           main = eval(expression(paste("tau^2",j))))
}

# --------------------------------------------
# Delta
# --------------------------------------------

par(mfrow = c(1,G), mar = c(5.1, 4.1, 4.1, 2.1))
for(j in 1:G){
 traceplot(mcmc(Delta.samp[,j]), xlab = "Iterações",
           main = eval(expression(paste("Delta",j))))
}

# --------------------------------------------
# nu
# --------------------------------------------

par(mfrow = c(1,1), mar = c(5.1, 4.1, 4.1, 2.1))
traceplot(mcmc(nu.samp), xlab = "Iterações",
          main = "nu")

# --------------------------------------------
# z
# --------------------------------------------

plot(z.samp[,1])
