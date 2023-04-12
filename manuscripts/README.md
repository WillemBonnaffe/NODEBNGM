# Bayesian neural gradient matching


## Overview

This repository contains a manuscript that introduces a novel method, Bayesian neural gradient matching (BNGM), to improve the speed and accuracy of neural ordinary differential equations (NODEs) fitting.

NODEs can be used to non-parametrically learn dynamical coupling between variables from time series data.
NODEs are slow to train, because we need to numerically solve the ODE system at every step of the optimisation.
BNGM, by interpolating the time series with neural networks, bypasses the numerical integration step.


## Method 

The method works in two steps (see figure below): 
* interpolating the time series data with neural networks (observation models), which provides an approximation of the dynamics and state of the variables 
* approximating the interpolated dynamics of the variables with neural networks (process model) which take as input the interpolated states

<p align="center">
<img align="middle" src="https://github.com/WillemBonnaffe/NODEBNGM/blob/main/examples/MEE_2023/.readme_files/fig_graphical_abstract_1.png" alt="" width="750" height="500" />
</p>


## Outcome

The outputs of the approach are:
* interpolations of the state and dynamics of the variables
* effects and contributions of each variables to the dynamics of the other variables

<p align="center">
<img align="middle" src="https://github.com/WillemBonnaffe/NODEBNGM/blob/main/examples/MEE_2023/.readme_files/fig_graphical_abstract_2.png" alt="" width="750" height="500" />
</p>


## Repository structure

The repository contains all the files for the manuscript: 
* `manuscript/` folder with the latex files to compile the manuscript, figures, and references
* `manuscript/figures/` folder with the latex files to compile the figures
* `manuscript/figures/scripts` with the R scripts that generate the results used for the figures of the manuscript

<p align="center">
<img align="middle" src="https://github.com/WillemBonnaffe/NODEBNGM/blob/main/examples/MEE_2023/.readme_files/fig_repo_overview.png" alt="" width="750" height="250" />
</p>


## Installation

`Latex` is required to compile the latex files of the manuscript and figures and `R` is required to run the scripts.


## Running the scripts

Instructions on how to run the scripts that generate the results used for the figures of the manuscript can be found in the `README.md` file of the `manuscript/figures/scripts` folder.
