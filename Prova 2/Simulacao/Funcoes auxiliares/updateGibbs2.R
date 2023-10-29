# deve ser rodado apos updateGibbs2
# atualiza somente parametros que dependem da componente
# param: lista com valores na ultima atualizacao dos parametros
#        deve estar nomeada na ordem c("prob","beta","tau2","Delta","alpha",
#        "nu","u","t","z")
# arg: tipo do parametro
# comp: qual componente atualizar
# lembrando que beta matriz pXG, tau2 e Delta vetores de tamanho G

updateGibbs2 = function(param,arg,comp){
 paramnames = names(param)
 if(all(arg != paramnames[2:4])){return(param)}
 novo = switch(arg,
               beta = full_beta(comp,param[["tau2"]],param[["Delta"]],param[["u"]],param[["t"]],param[["z"]]),
               tau2 = full_tau2(comp,param[["beta"]],param[["Delta"]],param[["u"]],param[["t"]],param[["z"]]),
               Delta = full_Delta(comp,param[["beta"]],param[["tau2"]],param[["u"]],param[["t"]],param[["z"]]))
 if(arg == "beta"){
  param[[arg]][,comp] = novo
 }else{
  param[[arg]][comp] = novo
 }
 return(param) 
}
