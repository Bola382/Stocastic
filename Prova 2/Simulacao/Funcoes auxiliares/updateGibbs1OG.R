# atualiza somente parametros que nao dependem da componente com a priori original para nu
# param: lista com valores na ultima atualizacao dos parametros
#        deve estar nomeada na ordem c("prob","beta","tau2","Delta","alpha",
#        "nu","u","t","z")
# arg: tipo do parametro
updateGibbs1OG = function(param,arg){
 paramnames = names(param)
 if(any(arg == paramnames[2:4])){return(param)}
 novo = switch(arg,
               prob = full_prob(param[["z"]]),
               alpha = full_alpha2(param[["nu"]]),
               nu = full_nu2(param[["prob"]],param[["beta"]],param[["tau2"]],param[["Delta"]],param[["alpha"]],param[["nu"]]),
               u = full_U(param[["beta"]],param[["tau2"]],param[["Delta"]],param[["nu"]],param[["t"]],param[["z"]]),
               t = full_T(param[["beta"]],param[["tau2"]],param[["Delta"]],param[["u"]],param[["z"]]),
               z = full_Z(param[["prob"]],param[["beta"]],param[["tau2"]],param[["Delta"]],param[["nu"]]))
 param[[arg]] = novo
 return(param)
}