# Loading packages
rm(list=ls())

source('funcoes_aula.R')

library(dplyr)
library(sommer)
library(AGHmatrix)
library(pheatmap)
library(EnvRtype)
library(ggplot2)

################################################################################
######################### Loading the data #####################################
################################################################################

####
## Phenotypic data - maize yield (tons per ha)

# Adjusted means for 100 genotypes across 5 environments
data <- read.csv('pheno_data.csv', sep=';', stringsAsFactors = T)
str(data)

####
## Genotypic data
load('Markers.Rda')

# 34,326 SNPs for 570 individuals
dim(M)

# Filtering data
gids <- levels(data$gid)

M <- M[gids,]
dim(M)

# Inspecting the matrix
M[1:10,1:10]

################################################################################
######################  Single-environment Genomic Selection ###################
################################################################################

# Let's filter environment "1_AN_LN"
data1 <- data[data$env == "1_AN_LN",] %>% droplevels()
str(data1)

# Building the G matrix -> VanRaden method
G <- Gmatrix(M)
dim(G)

# Heatmap
pheatmap(G)

####
## Building the model

# y = u + g + e --> means already adjusted

mmes.G1 <- mmes(value ~ 1,                  # y = [...]
                random = ~ vsm(ism(gid), Gu = G),  # [...] + g + e
                data = data1, getPEV = T)

summary(mmes.G1)

####
## Dissecting the model

# Intercept (u)
mmes.G1$b

# Genomic values (g)
(blups <- mmes.G1$u)

# Residuals (e)
mmes.G1$residuals

# Predictive Error Variance (PEV)
(pevs <- as.numeric(mmes.G1$uPevList[[1]]))  # In the same order as blups

####
## Visualizing genotypes

df.blups <- data.frame(blups, pevs, gid = rownames(blups))

ggplot(df.blups, aes(x = blups, y = reorder(as.factor(gid), blups))) +
  geom_point() +
  geom_errorbar(aes(xmax = blups + sqrt(pevs)*2, xmin = blups - sqrt(pevs)*2)) +
  labs(x = 'Genomic Value', y = 'gid')

####
## Selecting top performers

# Sorting
df.blups <- df.blups[order(df.blups$blups, decreasing = T),]

# Selecting the top 10
df.blups[1:10,]
mean(df.blups[1:10,]$blups)  # Genetic gain

################################################################################
######################  Genomic Prediction Single Env ##########################
################################################################################

# Can the model predict new genotypes or enable early selection?
# Suggestion: 5-fold cross-validation

set.seed(10)
gen.fold <- fold(5, 100)

(acc <- mmes.fold(dados = data1, fold = gen.fold, G = G))

mean(acc)
# Accuracy ~0.30

################################################################################
######################  Multi-environment Genomic Selection ####################
################################################################################

# Minimum: 4 environments

####
## Building the model

# y = u + a + g + gxa + e --> replicates already adjusted

mmes.G2 <- mmes(value ~ 1,  
                random = ~ env + vsm(ism(gid), Gu = G) + env:gid, 
                data = data, getPEV = T)

summary(mmes.G2)

####
## Dissecting the model

# Intercept (u)
mmes.G2$b 

# Genomic values (g) and G×E
mmes.G2$u

(blups <- mmes.G2$uList[2]$`vsm(ism(gid), Gu = G)`)

# Residuals (e)
mmes.G2$residuals

# Predictive Error Variance (PEV)
(pevs <- as.numeric(mmes.G2$uPevList[[1]]))

####
## Visualizing genotypes

df.blups <- data.frame(blups = blups, pevs, gid = rownames(blups))

ggplot(df.blups, aes(x = blups, y = reorder(as.factor(gid), blups))) +
  geom_point() +
  geom_errorbar(aes(xmax = blups + sqrt(pevs)*2, xmin = blups - sqrt(pevs)*2)) +
  labs(x = 'Genomic Value', y = 'Hybrid')

####
## Selecting top performers

# Sorting
df.blups <- df.blups[order(df.blups$mu, decreasing = T),]

# Top 10 selection
df.blups[1:10,]
mean(df.blups[1:10,]$mu)  # Genetic gain

################################################################################
######################  Multi-environment Genomic Prediction ###################
################################################################################

# Accuracy: correlation(y, yHat)

# CV2 - Random data partition
acc.cv2 <- mmesCV(dados = data, G = G, CV = '2', fold.n = 5)
colMeans(acc.cv2)  # ~0.56

# CV1 - Predicting genotypes
acc.cv1 <- mmesCV(dados = data, G = G, CV = '1', fold.n = 5)
colMeans(acc.cv1)  # ~0.50

# CV0 - Predicting environments
acc.cv0 <- mmesCV(dados = data, G = G, CV = '0', fold.n = 5)
colMeans(acc.cv0)  # ~0.43

# CV00 - Predicting genotypes and environments
acc.cv00 <- mmesCV(dados = data, G = G, CV = '00', fold.n = 5)
colMeans(acc.cv00)  # ~0.2

#####
## Adding covariance structure for G×E interaction

# CV2
acc.cv2 <- mmesCV(dados = data, G = G, CV = '2', fold.n = 5, covGE = T)
colMeans(acc.cv2)  # ~0.57

# CV1
acc.cv1 <- mmesCV(dados = data, G = G, CV = '1', fold.n = 5, covGE = T)
colMeans(acc.cv1)  # ~0.52

# CV0
acc.cv0 <- mmesCV(dados = data, G = G, CV = '0', fold.n = 5, covGE = T)
colMeans(acc.cv0)  # ~0.40

# CV00
acc.cv00 <- mmesCV(dados = data, G = G, CV = '00', fold.n = 5, covGE = T)
colMeans(acc.cv00)  # ~0.18

###########################################################
### Including environmental covariance structure ##########
###########################################################

# Importing environmental data
E <- readRDS('W_matrix_USP_248')
E[1:8, 1:4]
dim(E)

#Note: The covariate matrix is commonly referred to as W, and the similarity matrix as Omega.
#      Here, we chose to use E and W, respectively.

# Filtering environments
envs <- levels(data$env)
E <- E[envs,]

# Building the environmental similarity matrix
W <- gaussian(E)

pheatmap(W)  # Heatmap

# Genomic prediction

# CV2
acc.cv2 <- mmesCV(dados = data, G = G, CV = '2', fold.n = 5, covGE = T, W = W)
colMeans(acc.cv2)  # ~0.60

# CV1
acc.cv1 <- mmesCV(dados = data, G = G, CV = '1', fold.n = 5, covGE = T, W = W)
colMeans(acc.cv1)  # ~0.55

# CV0
acc.cv0 <- mmesCV(dados = data, G = G, CV = '0', fold.n = 5, covGE = T, W = W)
colMeans(acc.cv0)  # ~0.42

# CV00
acc.cv00 <- mmesCV(dados = data, G = G, CV = '00', fold.n = 5, covGE = T, W = W)
colMeans(acc.cv00)  # ~0.23

################################################################################
##########################  Obtaining environmental data #######################
################################################################################
# First, geographic coordinates and planting/harvest dates are required
# Simulated coordinates and dates for different locations in Brazil

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


# Using NASA POWER

# Loop over environments
for (i in seq_along(env)) {
  if (i == 1) df <- data.frame()
  
  aux <- get_power(
    community = "ag",
    lonlat = c(lon[i], lat[i]),
    pars = c("T2M_MAX", "T2M", "PRECTOTCORR", "T2M_MIN", "RH2M", "WS2M", 
             "CLRSKY_SFC_SW_DWN", "ALLSKY_SFC_SW_DWN"),
    dates = c(plant.date[i], harv.date[i]),
    temporal_api = "daily"
  )
  
  aux$daysFromStart <- 1:nrow(aux)
  aux$env <- env[i]
  
  df <- rbind(aux, df)
}

# Option: Processing by phenological phase (EnvRtype package)
id.var <- names(df)[c(8:16)]  # variable names
E <- W_matrix(env.data = data.frame(df), var.id = id.var,
              statistic = 'quantile', by.interval = TRUE, 
              time.window = c(0, 14, 35, 65, 90, 120))
E

# Using EnvRtype
df <- get_weather(lat = lat, lon = lon, 
                  start.day = plant.date, 
                  end.day = harv.date, 
                  env.id = env)

id.var <- names(df)[c(9:17)]  # variable names
E <- W_matrix(env.data = data.frame(df), var.id = id.var,
              statistic = 'quantile', by.interval = TRUE,
              time.window = c(0, 14, 35, 65, 90, 120))
E

# Similarity matrix
W <- gaussian(E)
pheatmap::pheatmap(W)
