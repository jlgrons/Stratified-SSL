# Stratified-SSL: Semi-supervised Learning under Stratified Sampling


## Note!

**Please use the new package [here](https://github.com/jlgrons/stratifiedSSL).**

## Overview 

R package and simulation studies for [Efficient Estimation and Evaluation of Prediction Rules in Semi-Supervised Settings under Stratified Sampling](https://arxiv.org/abs/2010.09443) by Jessica Gronsbell, Molei Liu, Lu Tian, and Tianxi Cai. This repo contains the following subfolders:

* stratifiedSSL: R functions for the package implementing the methods proposed in the manuscript.

* Simulation Studies: R code and directions to replicate the simulation studies in the paper. 
  * run_simulation.R: Primary file to run the simulations in the main text + supplement.
  * run_IE_simulation.R: File to run the simulation for the intrinsic efficient semi-supervised estimator.
  * IE_helper_functions.R: File with helper functions for the intrinisic efficient simulations.
  
* Other Files: Previous versions of the R code for debugging purposes and documentation.

Code is maintained by Jessica Gronsbell and Molei Liu.

## Installation

```r
# install.packages("remotes")
remotes::install_github("jlgrons/Stratified-SSL/stratifiedSSL")
```
Required packages: MASS, evd, glmnet, matrixStats, methods, randomForest, stepPlr.

## Getting Started

Install the package as above and then run an example with 'run_simulation.R' in the Simulation Studies folder.
