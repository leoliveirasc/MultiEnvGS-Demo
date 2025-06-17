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

