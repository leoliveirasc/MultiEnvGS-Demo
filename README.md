# MultiEnvGS-Demo
This repository contains practical R scripts to demonstrate genomic selection models across multiple environments, with a focus on applications in plant breeding education and training.

## What's Included (English)
1) Cross-validation functions: Scripts implementing different cross-validation schemes (CV0, CV00, CV1, and CV2) commonly used in genomic selection.
2) Bilingual instructions: A script with usage instructions in both Portuguese (pt-BR) and English (ENG).
3) Example datasets: Phenotypic, genomic, and environmental data provided via https://github.com/gcostaneto/KernelMethods.
4) Modeling code: Scripts using the sommer package to perform multi-environment genomic selection.

**Model specification**: The models used follow the structure:  
`y = u + e + g + g×e`  

Where:  
- `y`: phenotype  
- `u`: overall mean  
- `e`: environmental effect or environmental covariates effects  
- `g`: genotypic or additive genomic effect  
- `g×e`: genotype-by-environment interaction 


##  O que está incluído (Português)
1) Funções de validação cruzada: Scripts que implementam diferentes esquemas de validação cruzada (CV0, CV00, CV1 e CV2), comumente usados na seleção genômica.
2) Instruções bilíngues: Script com orientações de uso em Português (pt-BR) e Inglês (ENG).
3) Conjuntos de dados de exemplo: Dados fenotípicos, genômicos e ambientais disponíveis em https://github.com/gcostaneto/KernelMethods
4) Código para modelagem: Scripts que utilizam o pacote sommer para análise de seleção genômica em múltiplos ambientes.

**Especificação do modelo**: Os modelos seguem a estrutura:
y = u + e + g + g×e

Onde:
- `y`: fenótipo
- `u`: média geral
- `e`: efeito do ambiente ou efeitos de covariáveis ambientais
- `g`: efeito genotípico ou efeito genômico aditivo
- `g×e`: interação genótipo × ambiente



##  Suggested References / Referências Sugeridas
- Jarquín, D., Crossa, J., Lacaze, X., Du Cheyron, P., Daucourt, J., Lorgeou, J., Piraux, F., Guerreiro, L., Pérez, P., Calus, M., Burgueño, J., de los Campos, G., 2014. A reaction norm model for genomic selection using high-dimensional genomic and environmental data. Theoretical and Applied Genetics 127, 595–607. https://doi.org/10.1007/s00122-013-2243-1

- Costa-Neto, G., Fritsche-Neto, R., Crossa, J., 2021a. Nonlinear kernels, dominance, and envirotyping data increase the accuracy of genome-based prediction in multi-environment trials. Heredity (Edinb) 126, 92–106. https://doi.org/10.1038/s41437-020-00353-1

- Covarrubias-Pazaran, G, 2025. Quantitative genetics using the sommer package. Available at: https://cran.r-project.org/web/packages/sommer/vignettes/sommer.qg.html
