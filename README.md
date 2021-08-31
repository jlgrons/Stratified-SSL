# Stratified-SSL: Semi-supervised learning under stratified sampling

R package and simulation studies for [Efficient Estimation and Evaluation of Prediction Rules in Semi-Supervised Settings under Stratified Sampling](https://arxiv.org/abs/2010.09443) by Jessica Gronsbell, Molei Liu, Lu Tian, and Tianxi Cai.

This repo contains the following subfolders:

* stratifiedSSL: R functions for the package implementing the methods proposed in the manuscript.
* Simulation Studies: R code and directions to replicate the simulation studies in the paper. 
* oldCode: Old versions of the R code for debugging purposes.

## Installation

```r
# install.packages("remotes")
remotes::install_github("jlgrons/Stratified-SSL/stratifiedSSL")
```

## Getting Started

You can run an example in 'run_simulation.R' in the Simulation Studies folder.
