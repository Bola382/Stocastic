# ==================================================================
# Verifica a convergência de uma cadeia MCMC
# as cadeias geradas nao estao disponiveis por conta de seu tamanho
# ==================================================================
library(coda)

GplusHat = ncol(prob_ok)

# --------------------------------------------
# prob
# --------------------------------------------

par(mfrow = c(1, GplusHat), mar = c(5.1, 4.1, 4.1, 2.1))
for (j in 1:GplusHat) {
 traceplot(mcmc(prob_ok[, j]), xlab = "Iterações",
           main = eval(expression(paste("p", j))))
}
for (j in 1:GplusHat) {
 acf(prob_ok[, j], xlab = "Iterações",
     main = eval(expression(paste("p", j))))
}

# --------------------------------------------
# Beta
# --------------------------------------------
p = ncol(X)

par(mfrow = c(p, GplusHat), mar = c(4, 4, 3, 1))
for (i in 1:p) {
 for (j in 1:GplusHat) {
  traceplot(mcmc(beta_ok[, i, j]), xlab = "Iterações",
            main = eval(expression(paste("beta", i, j))))
 }
}

# --------------------------------------------
# tau2
# --------------------------------------------

par(mfrow = c(1, GplusHat), mar = c(5.1, 4.1, 4.1, 2.1))
for (j in 1:GplusHat) {
 traceplot(mcmc(tau2_ok[, j]), xlab = "Iterações",
           main = eval(expression(paste("tau^2", j))))
}

for (j in 1:GplusHat) {
 acf(tau2_ok[, j], xlab = "Iterações",
     main = eval(expression(paste("tau^2", j))))
}

# --------------------------------------------
# Delta
# --------------------------------------------

par(mfrow = c(1, GplusHat), mar = c(5.1, 4.1, 4.1, 2.1))
for (j in 1:GplusHat) {
 traceplot(mcmc(Delta_ok[, j]), xlab = "Iterações",
           main = eval(expression(paste("Delta", j))))
}

for (j in 1:GplusHat) {
 acf(Delta_ok[, j], xlab = "Iterações",
     main = eval(expression(paste("Delta", j))))
}

# --------------------------------------------
# nu
# --------------------------------------------

par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))
traceplot(mcmc(nu_ok), xlab = "Iterações",
          main = "nu")

acf(nu_ok, xlab = "Iterações", main = "nu")

# --------------------------------------------
# z
# --------------------------------------------

plot(z_ok[, 1])
