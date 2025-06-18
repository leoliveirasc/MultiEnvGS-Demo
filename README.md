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
#' Randomly split observations into (approximately) equal folds
#'
#' @description
#' Utility to generate indices for \(k\)-fold cross‑validation.  
#' The function shuffles the integers \eqn{1:n} and returns a list whose
#' elements are integer vectors with the positions assigned to each fold.
#' Folds differ by at most one observation.
#'
#' @param fold.n Integer ≥ 2. Desired number of folds.
#' @param n      Integer ≥ fold.n. Total number of observations to split.
#'
#' @return A list of length \code{fold.n}; each element is a vector of indices.
#'
#' @examples
#' set.seed(42)
#' folds <- fold(5, 103)
#' lengths(folds)             # numbers per fold
#' intersect(folds[[1]], folds[[2]])  # empty — no overlap
#'
#' @export
fold <- function(fold.n, n) { ... }

##  Suggested References / Referências Sugeridas
- Jarquín, D., Crossa, J., Lacaze, X., Du Cheyron, P., Daucourt, J., Lorgeou, J., Piraux, F., Guerreiro, L., Pérez, P., Calus, M., Burgueño, J., de los Campos, G., 2014. A reaction norm model for genomic selection using high-dimensional genomic and environmental data. Theoretical and Applied Genetics 127, 595–607. https://doi.org/10.1007/s00122-013-2243-1

- Costa-Neto, G., Fritsche-Neto, R., Crossa, J., 2021a. Nonlinear kernels, dominance, and envirotyping data increase the accuracy of genome-based prediction in multi-environment trials. Heredity (Edinb) 126, 92–106. https://doi.org/10.1038/s41437-020-00353-1

- Covarrubias-Pazaran, G, 2025. Quantitative genetics using the sommer package. Available at: https://cran.r-project.org/web/packages/sommer/vignettes/sommer.qg.html
