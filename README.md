# MultiEnvGS-Demo
This repository contains practical R scripts to demonstrate genomic selection models across multiple environments, with a focus on applications in plant breeding education and training.

## What's Included
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


## Suggested references
- Jarquín, D., Crossa, J., Lacaze, X., Du Cheyron, P., Daucourt, J., Lorgeou, J., Piraux, F., Guerreiro, L., Pérez, P., Calus, M., Burgueño, J., de los Campos, G., 2014. A reaction norm model for genomic selection using high-dimensional genomic and environmental data. Theoretical and Applied Genetics 127, 595–607. https://doi.org/10.1007/s00122-013-2243-1

- Costa-Neto, G., Fritsche-Neto, R., Crossa, J., 2021a. Nonlinear kernels, dominance, and envirotyping data increase the accuracy of genome-based prediction in multi-environment trials. Heredity (Edinb) 126, 92–106. https://doi.org/10.1038/s41437-020-00353-1

- Covarrubias-Pazaran, G, 2025. Quantitative genetics using the sommer package. Available at: https://cran.r-project.org/web/packages/sommer/vignettes/sommer.qg.html
