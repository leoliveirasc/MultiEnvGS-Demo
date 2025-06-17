# Lista de pacotes CRAN
cran_pkgs <- c("dplyr", "sommer", "AGHmatrix", "pheatmap", "ggplot2")

# Instala os pacotes do CRAN se necessário
for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

# Instala e carrega EnvRtype (GitHub)
if (!requireNamespace("EnvRtype", quietly = TRUE)) {
  if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")
  }
  library(devtools)
  install_github('allogamous/EnvRtype', force = TRUE)
}
library(EnvRtype)




fold<-function(fold.n, n){
  
  r1<-round(n/fold.n)
  r2<- n- r1*(fold.n-1)
  
  return(split(sample(1:n,n),rep(1:fold.n,c(rep(r1,fold.n-1),r2) )))
  
  
}


mmes.fold<-function(dados, G,fold){
  
  # fold<-gen.fold
  # dados<-data1

  
  for(i in 1:length(fold)){
    if(i==1) acc<-vector()
    
    data.fold<-dados
    tst<-fold[[i]]
    
    data.fold$value[tst]<-NA
    mmes.fold <- mmes(value ~ 1,   #y
                      random = ~ vsm(ism(gid), Gu = G), # g + e
                      data = data.fold,getPEV = T)
    
    
   acc<-c(acc, cor(mmes.fold$u[tst],dados$value[tst]))
    
  }
  return(acc)
}


mmesCV<-function(dados, G = NULL, W = NULL,CV=NULL,fold.n=5, looEnv=T,covGE=F){
  set.seed(100)

  # fold.n<-5; looEnv=T;covGE=T
  # dados<-data
  # CV<-'0'
  # G=G;W = E

  ########
  ## Cross Validation Schemes
  
  if(CV=='2') fold.i<-fold(fold.n, nrow(dados))
  if(CV=='1'){
    gids<-levels(dados$gid)
    gen.fold<-fold(fold.n, length(gids))
    
    for(i in 1:length(gen.fold)){
      if(i==1)fold.i<-list()
      fold.i[[i]]<-which(dados$gid %in% gids[gen.fold[[i]]])
    }
    
    
  }

  if(CV=='0'){
    envs<-levels(dados$env)
    
    if(looEnv==T) fold.env=length(envs) else fold.env=fold.n
    env.fold<-fold(fold.env, length(envs))
    
    for(i in 1:length(env.fold)){
      if(i==1)fold.i<-list()
      fold.i[[i]]<-which(dados$env %in% envs[env.fold[[i]]])
    }
  }
  
  if(CV=='00'){
    envs<-levels(dados$env)
    
    if(looEnv==T) fold.env=length(envs) else fold.env=fold.n
    env.fold<-fold(fold.env, length(envs))
    
    for(i in 1:length(env.fold)){
      env.fold[[i]]<-which(dados$env %in% envs[env.fold[[i]]])
    }
    
    gids<-levels(dados$gid)
    gen.fold<-fold(fold.n, length(gids))
    
    for(i in 1:length(gen.fold)){
      gen.fold[[i]]<-which(dados$gid %in% gids[gen.fold[[i]]])
    }
    
    
    df.cv00<-expand.grid(env=1:fold.env,gid=1:fold.n)
    
    
    for(i in 1:nrow(df.cv00)){
      if(i==1)fold.i<-list()
      fold.i[[i]]<-union(gen.fold[[df.cv00$gid[i]]],env.fold[[df.cv00$env[i]]])
      
    }
    
  }
  ########
  # Kinship Matrices
  
  if(is.null(G)==T){
    gids<-levels(dados$gid)
    G<-diag(length(gids))
    rownames(G)<-colnames(G)<-gids
    
  }

  if(is.null(W)==T){
    envs<-levels(dados$env)
    W<-diag(length(envs))
    rownames(W)<-colnames(W)<-envs
    
  }
  
  if(covGE==F){
    
    gids<-levels(dados$gid)
    Ig<-diag(length(gids))
    
    envs<-levels(dados$env)
    Iw<-diag(length(envs))
    
    df.comp<-expand.grid(env = unique(dados$env), gid = unique(dados$gid),value=1)
    dim(df.comp)
    
    Ze<-model.matrix(value ~ -1 + env, data=df.comp)
    Zg<-model.matrix(value ~ -1 + gid, data=df.comp)
    
    ZEZ<-Ze%*%Iw%*%t(Ze)
    ZGZ<-Zg%*%Ig%*%t(Zg)
    
    GE<-ZEZ*ZGZ
    rownames(GE)<-colnames(GE)<-paste(dados$env,dados$gid,sep=':')
    
  }else{
    
    df.comp<-expand.grid(env = unique(dados$env), gid = unique(dados$gid),value=1)
    dim(df.comp)
    
    Ze<-model.matrix(value ~ -1 + env, data=df.comp)
    Zg<-model.matrix(value ~ -1 + gid, data=df.comp)
    
    ZEZ<-Ze%*%W%*%t(Ze)
    ZGZ<-Zg%*%G%*%t(Zg)
    
    GE<-ZEZ*ZGZ
    rownames(GE)<-colnames(GE)<-paste(dados$env,dados$gid,sep=':')

  }

  
  ########
  # Validation
  
  for(i in 1:length(fold.i)){
    if(i==1) acc<-acc.gge<-vector()
    
    data.fold<-dados
    tst<-fold.i[[i]]
    
    data.fold$value[tst]<-NA
    mmes.fold <- mmes(value ~ 1,   #y = u +[...]
                      random = ~ vsm(ism(env), Gu = W) # [...] + a
                      + vsm(ism(gid), Gu = G)# [...] + g
                      + vsm(ism(env:gid), Gu = GE), # [...] + gxa
                      data = data.fold)
    
    
    # Intercepto
    mu <- mmes.fold$b[1, 1]  # valor numérico
    
    # Efeitos aleatórios
    blup_env <- mmes.fold$uList[[1]]
    blup_gid <- mmes.fold$uList[[2]]
    blup_ge  <- mmes.fold$uList[[3]]
    
    #Renaming
    rownames(blup_gid)<-sort(gids)
    rownames(blup_ge)<-sort(paste(dados$env,dados$gid,sep=':'))

    #yHat
    data.fold$yhat <- mu +
      blup_env[data.fold$env, 1] +
      blup_gid[data.fold$gid, 1] +
      blup_ge[paste(data.fold$env, data.fold$gid, sep = ":"), 1]
    
    #gge
    data.fold$gge <-
      blup_gid[data.fold$gid, 1] +
      blup_ge[paste(data.fold$env, data.fold$gid, sep = ":"), 1]
    
  if(CV=="00"){
    
    tst.cv00<-intersect(gen.fold[[df.cv00$gid[i]]],env.fold[[df.cv00$env[i]]])
    
    acc<-c(acc, cor(data.fold$yhat[tst.cv00],dados$value[tst.cv00]))
    acc.gge<-c(acc.gge, cor(data.fold$gge[tst.cv00],dados$value[tst.cv00]))
    
  }else{
    
    acc<-c(acc, cor(data.fold$yhat[tst],dados$value[tst]))
    acc.gge<-c(acc.gge, cor(data.fold$gge[tst],dados$value[tst]))
    
  }

    
  }
  return(data.frame(acc,acc.gge))
}

