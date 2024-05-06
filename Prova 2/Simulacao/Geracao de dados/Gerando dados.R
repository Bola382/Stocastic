#=======================================================================
# Gerando um conjunto de dados do modelo de misturas de distribuicao St
#=======================================================================

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list=ls())

library("rgl")
source("rregmixst.R")

tmp_n = 1000
probs = c(.8,.15,.05)
beta = cbind(c(0,4), c(-10,-5), c(-25,-5))
sigma = c(1,2,3)
lambda = c(10,-8,-8)
nu = 3

set.seed(16)

tmp_x1 = rnorm(tmp_n)

tmp_covs = cbind(1, tmp_x1)

tmp_resul = rregmixst(tmp_covs, probs, beta, sigma, lambda, nu)

data = data.frame(group = tmp_resul$grupo, resp = tmp_resul$resposta)
data = cbind.data.frame(data,tmp_covs[,-1])
colnames(data) = c("comp", "resp", "x1")

plot(data$x1, data$resp, col = data$comp, cex = .9, pch = 16)

rm(list = ls(pattern="tmp_"))
rm(list = lsf.str())

save.image("n1000type3.RData")
