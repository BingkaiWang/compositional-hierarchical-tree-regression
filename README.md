# R code for Regularized regression on compositional trees with application to MRI analysis

## Overview

This repo stores R code for simulations and data analysis in the [paper](https://arxiv.org/abs/2104.07113).

## Data analysis

-  __data-analysis/ADNI-analysis.R__ stores the code for the MRI anlaysis. The MRI data is available upon request from from the Alzheimer's Disease Neu-
roimaging Initiative (ADNI) database (adni.loni.usc.edu).

## Simulation

- __scenario1.R - simulation4.R__ correspond to Scenarios 1-4 of Section 5 of the paper.
- For Scenarios 3-4, the MRI data is stored at __sim2-data.rds__, which stores the compositional tree structure and volumetric data. This file can be obtained by running __data-analysis/ADNI-analysis.R__ given that the MRI data has been acquired.
