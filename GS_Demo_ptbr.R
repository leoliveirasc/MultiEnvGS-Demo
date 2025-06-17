#Carregando pacotes
rm(list=ls())

source('funcoes_aula.R')

library(dplyr)
library(sommer)
library(AGHmatrix)
library(pheatmap)
library(EnvRtype)
library(ggplot2)

################################################################################
######################### Carregando os dados ##################################
################################################################################

####
## Dados fenot?picos - produtividade milho (tons per ha)

# Medias ajustadas de 100 genotipos para 5 ambientes
data<-read.csv('pheno_data.csv',sep=';',stringsAsFactors = T)
str(data)

####
## Dados genot?picos
load('Markers.Rda')

#34326 SNPs para 570 individuos
dim(M)

#filtrando dados
gids<-levels(data$gid)

M<-M[gids,]
dim(M)

#Investigando a matriz
M[1:10,1:10]

################################################################################
######################  Sele??o Gen?mica 1 ambiente ############################
################################################################################

#Vamos filtrar o ambiente  "1_AN_LN"
data1<-data[data$env=="1_AN_LN",] %>% droplevels()
str(data1)

#Construindo a matriz G -> VanRaden

G<-Gmatrix(M)
dim(G)

#Heatmap
pheatmap(G)

####
## Construindo o modelo

# y = u + g + e --> M?dias j? diluidas

mmes.G1 <- mmes(value ~ 1,   #y = [...]
               random = ~ vsm(ism(gid), Gu = G), # [...] + g + e
               data = data1,getPEV = T)

summary(mmes.G1)

####
## Dissecando o modelo

# Intercepto (u)
mmes.G1$b

#Valores reprodutivos (g)
(blups<- mmes.G1$u)

#res?duos (e)
mmes.G1$residuals

#Variancia do Erro Preditivo (PEV)
(pevs<-as.numeric(mmes.G1$uPevList[[1]])) # Na ordem dos blups

####
## Visualizando os gen?tipos

df.blups<- data.frame(blups,pevs,gid=rownames(blups))

ggplot(df.blups,aes(x=blups,y=reorder(as.factor(gid), blups)))+
  geom_point()+
  geom_errorbar(aes(xmax= blups + sqrt(pevs)*2,xmin= blups - sqrt(pevs)*2 )) +
  labs(x='Valor reprodutivo',y='gid')

####
## Selecionando os melhores

# Ordenando
df.blups<-df.blups[order(df.blups$blups,decreasing = T),]

# Selecionando os 10 melhores
df.blups[1:10,]
mean(df.blups[1:10,]$blups) # Ganho genetico


################################################################################
######################  Predicao Genomica 1 ambiente ###########################
################################################################################

#O modelo consegue fazer predicoes de novos genotipos ou selecao precoce?
#Sugestao validacao: 5 fold 

set.seed(10)
gen.fold<-fold(5,100)

(acc<-mmes.fold(dados = data1,fold = gen.fold,G = G))

mean(acc)
#acuracia ~0.30

################################################################################
######################  Selecao Genomica multi-ambiente ########################
################################################################################

#Minimo 4 ambientes!

####
## Construindo o modelo

# y = u + a + g + gxa + e --> repeticoes ja diluidas

mmes.G2 <- mmes(value ~ 1 ,   #y = u + [...]
                random = ~ env + vsm(ism(gid), Gu = G) + env:gid, # [...] a + g + gxa + e
                data = data,getPEV = T)


summary(mmes.G2)

####
## Dissecando o modelo

# Intercepto (u)
mmes.G2$b 

#Valores reprodutivos (g) e gxa 
mmes.G2$u

(blups<-mmes.G2$uList[2]$`vsm(ism(gid), Gu = G)`)

#res?duos (e)
mmes.G2$residuals

#Variancia do Erro Preditivo (PEV)
(pevs<-as.numeric(mmes.G2$uPevList[[1]])) # Na ordem dos blups

####
## Visualizando os gen?tipos

df.blups<- data.frame(blups=blups,pevs,gid=rownames(blups))

ggplot(df.blups,aes(x=blups,y=reorder(as.factor(gid), blups)))+
  geom_point()+
  geom_errorbar(aes(xmax= blups + sqrt(pevs)*2,xmin= blups - sqrt(pevs)*2 ))+
  labs(x='Valor reprodutivo',y='Hibrido')

####
## Selecionando os melhores

# Ordenando
df.blups<-df.blups[order(df.blups$mu,decreasing = T),]

# Selecionando os 10 melhores
df.blups[1:10,]
mean(df.blups[1:10,]$mu) # Ganho genetico

################################################################################
######################  Predicao Genomica multi-ambiente #######################
################################################################################

#Acuracia cor(y,yHat)

#CV2 - dados aleatorios
acc.cv2<-mmesCV(dados = data,G = G,CV = '2',fold.n = 5)
colMeans(acc.cv2) # ~ 0.56 

#CV1 - predicao de genotipos
acc.cv1<-mmesCV(dados = data,G = G,CV = '1',fold.n = 5)
colMeans(acc.cv1) # ~ 0.50 

#CV0 - predicao de ambientes
acc.cv0<-mmesCV(dados = data,G = G,CV = '0',fold.n = 5)
colMeans(acc.cv0) # ~ 0.43 

#CV00 - predicao de genotipos e ambientes
acc.cv00<-mmesCV(dados = data,G = G,CV = '00',fold.n = 5)
colMeans(acc.cv00) # ~ 0.2


#####
## Incluindo estruturas de covariancia para a interacao

#CV2
acc.cv2<-mmesCV(dados = data,G = G,CV = '2',fold.n = 5, covGE = T)
colMeans(acc.cv2) # ~ 0.57 

#CV1
acc.cv1<-mmesCV(dados = data,G = G,CV = '1',fold.n = 5, covGE = T)
colMeans(acc.cv1) # ~ 0.52 

#CV0
acc.cv0<-mmesCV(dados = data,G = G,CV = '0',fold.n = 5, covGE = T)
colMeans(acc.cv0) # ~ 0.40 

#CV00
acc.cv00<-mmesCV(dados = data,G = G,CV = '00',fold.n = 5, covGE = T)
colMeans(acc.cv00) # ~ 0.18 


###########################################################
### Incluindo estruturas de covariancia para o ambiente ###
###########################################################

#Importando dados ambientais
E<-readRDS('W_matrix_USP_248')
E[1:8,1:4]
dim(E)

#Obs: A matriz de covariaveis eh comumente chamada de W e a matriz de similaridade de Omega,
#     aqui optamos por E e W respectivamente.

#Filtrando ambientes
envs<-levels(data$env)
E<-E[envs,]

#Construindo a matriz de similaridade ambiental
W<-gaussian(E)

pheatmap(W) #Heatmap

#Predicao Genomica

#CV2
acc.cv2<-mmesCV(dados = data,G = G,CV = '2',fold.n = 5,covGE = T, W = W)
colMeans(acc.cv2) # ~ 0.60 

#CV1
acc.cv1<-mmesCV(dados = data,G = G,CV = '1',fold.n = 5,covGE = T, W = W)
colMeans(acc.cv1) # ~ 0.55

#CV0
acc.cv0<-mmesCV(dados = data,G = G,CV = '0',fold.n = 5,covGE = T, W = W)
colMeans(acc.cv0) # ~ 0.42 

#CV00
acc.cv00<-mmesCV(dados = data,G = G, CV = '00',fold.n = 5,covGE = T, W = W)
colMeans(acc.cv00) # ~ 0.23 


################################################################################
##########################  Obtendo dados ambientais ###########################
################################################################################

#Primeiro sao necessarias as coordenadas geograficas, data de inicio e final 
# Coordenadas e datas simuladas para diferentes locais no Brasil

lat <- c(-16.6869,   # Goiás (Goiânia)
         -12.6819,   # Mato Grosso (Sorriso)
         -15.6146,   # Mato Grosso (Rondonópolis)
         -22.9099,   # São Paulo (Campinas)
         -18.9096)   # Minas Gerais (Uberaba)

lon <- c(-49.2648,   # Goiás (Goiânia)
         -55.7200,   # Mato Grosso (Sorriso)
         -54.3411,   # Mato Grosso (Rondonópolis)
         -47.0626,   # São Paulo (Campinas)
         -47.0626)   # Minas Gerais (Uberaba)

env <- c("GO1", "MT1", "MT2", "SP1", "MG1")

plant.date <- c("2021-01-20",  # Goiás
                "2021-01-15",  # MT - Sorriso
                "2021-01-18",  # MT - Rondonópolis
                "2021-01-22",  # SP
                "2021-01-17")  # MG

harv.date <- c("2022-07-28",   # Goiás
               "2022-07-20",   # MT - Sorriso
               "2022-07-25",   # MT - Rondonópolis
               "2022-07-30",   # SP
               "2022-07-22")   # MG


#Com Nasa Power

# Loop por ambiente
for (i in seq_along(env)) {
  if(i==1) df<-data.frame()
  aux<- get_power(
    community = "ag",
    lonlat = c(lon[i], lat[i]),
    pars = c("T2M_MAX","T2M", "PRECTOTCORR","T2M_MIN", "RH2M", "WS2M","CLRSKY_SFC_SW_DWN", "ALLSKY_SFC_SW_DWN"),
    dates = c(plant.date[i], harv.date[i]),
    temporal_api = "daily"
  )
  
  aux$daysFromStart<-1:nrow(aux)
  aux$env<-env[i]
  

  df<-rbind(aux,df)
}

#Opcao: Processamento por fase fenológica
id.var <- names(df)[c(8:16)] # nome das vari?veis
E <-W_matrix(env.data = data.frame(df),var.id = id.var,
             statistic = 'quantile',by.interval = TRUE,time.window = c(0,14,35,65,90,120))
E

#Com EnvRtype
df<-get_weather(lat=lat,lon=lon,start.day =plant.date,end.day = harv.date,env.id = env )

id.var <- names(df)[c(9:17)] # nome das vari?veis
E <-W_matrix(env.data = data.frame(df),var.id = id.var,
             statistic = 'quantile',by.interval = TRUE,time.window = c(0,14,35,65,90,120))
E


#Matriz de similaridade
W<-gaussian(E)
pheatmap::pheatmap(W)

