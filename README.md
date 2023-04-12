# NODEBNGM
 

----------------
## Overview

This repository introduces a novel method, Bayesian neural gradient matching (BNGM), to improve the speed and accuracy of neural ordinary differential equations (NODEs) fitting.

NODEs can be used to nonparametrically learn dynamical coupling between variables from time series data. However, NODEs are slow to train because we need to numerically solve the ODE system at every step of the optimisation. BNGM, by interpolating the time series with neural networks, which allows us to fit NODEs directly to the interpolated dynamics, hence bypassing the numerical integration step. 

The two main advantage of NODEBNGM are:
* faster to fit than standard NODEs
* no parametric assumptions required to estimate interactions between variables


----------------
## Method 

### Principle 

The method works in two steps (see figure below): 
* interpolating the time series data with neural networks (observation models), which provides an approximation of the dynamics and state of the variables 
* approximating the interpolated dynamics of the variables with neural networks (process model) which take as input the interpolated states

<p align="center">
<img align="middle" src="https://github.com/WillemBonnaffe/NODEBNGM/blob/main/.readme_files/fig_graphical_abstract_1.png" alt="" width="750" height="500" />
</p>

### Inputs

The method uses time series datasets, which should have one time series per variable in the system.

<p align="center">
<img align="middle" src="https://github.com/WillemBonnaffe/NODEBNGM/blob/main/.readme_files/fig_time_series.png" alt="" width="500" height="500" />
</p>


### Outcomes

The outputs of the approach are:
* interpolations of the state and dynamics of the variables
* effects and contributions of each variables to the dynamics of the other variables

<p align="center">
<img align="middle" src="https://github.com/WillemBonnaffe/NODEBNGM/blob/main/.readme_files/fig_DIN.png" alt="" width="500" height="500" />
</p>


----------------
## Repository structure

The folder `manuscripts` contains all the files necessary for the manuscript and the preprint: 
* `preprint` folder with the latex files, figures, and scripts to compile the arXiv preprint
* `manuscript/` folder with the latex files to compile the manuscript, figures, and references
* `manuscript/figures/` folder with the latex files to compile the figures
* `manuscript/figures/scripts` with the R scripts that generate the results used for the figures of the manuscript

The folder includes an artificial time series of a three-species prey-predator system, 3 replicated time series of an algae, flagellate, and rotifer system, and the hare-lynx time series.

The folder `scripts` contains the template repository for fitting NODE by BNGM to new time series datasets.

Further instructions on how to run NODEs can be found in these directories.


----------------
## Installation

* `LaTeX` (TeX live 2022) is required to compile the latex files of the manuscript and figures 
* `R` (v4.0.0 or later) is required to run the scripts


----------------
## Running the scripts

Instructions on how to run the scripts that generate the results used for the figures of the manuscript can be found in the `README.md` file of the relevant script folders.
