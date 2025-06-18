# MultiEnvGS-Demo
This repository contains practical R scripts to demonstrate genomic selection models across multiple environments, with a focus on applications in plant breeding education and training.

## What's Included (English)
1) Cross-validation functions: Scripts implementing different cross-validation schemes (CV0, CV00, CV1, and CV2) commonly used in genomic selection.
2) Bilingual instructions: Scripts with usage instructions available in both Portuguese (GS_Demo_ptbr.R) and English (GS_Demo_eng.R), including modeling code that uses the sommer package to perform multi-environment genomic selection.
3) Example datasets: Phenotypic (pheno_data.csv), genomic (Markers.Rda), and environmental data (W_matrix_USP_248) provided via https://github.com/gcostaneto/KernelMethods.

**Model specification**: The models used follow the structure:  
`y = u + e + g + g×e`  

Where:  
- `y`: phenotype  
- `u`: intercept  
- `e`: environmental effect or environmental covariates effects  
- `g`: genotypic or additive genomic effect  
- `g×e`: genotype-by-environment interaction 


##  O que está incluído (Português)
1) Funções de validação cruzada: Scripts que implementam diferentes esquemas de validação cruzada (CV0, CV00, CV1 e CV2), comumente usados na seleção genômica.
2) Instruções bilíngues: Scripts com instruções de uso disponíveis em português (GS_Demo_ptbr.R) e inglês (GS_Demo_eng.R), incluindo código de modelagem que utiliza o pacote sommer para realizar seleção genômica multiambiente.
3) Conjuntos de dados de exemplo: Dados fenotípicos (pheno_data.csv), genômicos (Markers.Rda) e ambientais (W_matrix_USP_248) disponíveis em https://github.com/gcostaneto/KernelMethods

**Especificação do modelo**: Os modelos seguem a estrutura:
`y = u + e + g + g×e`  

Onde:
- `y`: fenótipo
- `u`: intercepto
- `e`: efeito do ambiente ou efeitos de covariáveis ambientais
- `g`: efeito genotípico ou efeito genômico aditivo
- `g×e`: interação genótipo × ambiente

## Pipeline
<img src="Pipeline.png" alt="Esquema geral da análise" width="700">

## Functions
### fold
fold(fold.n, n)

##### fold.n: number of partitions
##### n: number of observations

### mmes.fold
for single environments  
mmes.fold(dados, G, fold)

##### dados: phenotypic data containing 'gid' and 'env' columns  
##### G: genomic additive matrix  
##### fold: list obtained with fold() function  

### mmesCV
for multi-environments  
mmesCV(dados, G = NULL, W = NULL, CV = NULL, fold.n = 5, looEnv = TRUE, covGE = FALSE)

##### dados: phenotypic data containing 'gid' and 'env' columns  
##### G: genomic additive matrix  
##### W: environmental similarity matrix  
##### fold.n: number of partitions in CV. The default is 5.  
##### looEnv: for CV = '00' only; indicates if leave-one-out cross-validation should be applied specifically over environments, while genotypes are partitioned into folds defined by fold.n  
##### covGE: if covariances for GxE will be constructed (Jarquín et al., 2014)  


##  Suggested References / Referências Sugeridas
- Jarquín, D., Crossa, J., Lacaze, X., Du Cheyron, P., Daucourt, J., Lorgeou, J., Piraux, F., Guerreiro, L., Pérez, P., Calus, M., Burgueño, J., de los Campos, G., 2014. A reaction norm model for genomic selection using high-dimensional genomic and environmental data. Theoretical and Applied Genetics 127, 595–607. https://doi.org/10.1007/s00122-013-2243-1

- Costa-Neto, G., Fritsche-Neto, R., Crossa, J., 2021a. Nonlinear kernels, dominance, and envirotyping data increase the accuracy of genome-based prediction in multi-environment trials. Heredity (Edinb) 126, 92–106. https://doi.org/10.1038/s41437-020-00353-1

- Covarrubias-Pazaran, G, 2025. Quantitative genetics using the sommer package. Available at: https://cran.r-project.org/web/packages/sommer/vignettes/sommer.qg.html
