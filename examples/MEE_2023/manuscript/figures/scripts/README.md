# Scripts 


## Aim

This document provides an overview of the script folder and the directories that it contains.
This folder contains all the scripts and data necessary to replicate the results presented in the manuscript.


## Repository structure

Each folder contains the scripts of a specific analysis of the manuscript:
* `benchmark` contains the scripts to compare the runtime and accuracy of fitting NODEs by BNGM compared to standard NODEs, parametric polynomical ODEs, and convergent cross mapping (CCM).
* `cmdline_HL` contains the scripts to analyse of the hare-lynx time series dataset.
* `cmdline_TS1-3` contain the scripts to analyse the three replicated time series of the algae, flagellate, rotifer tri-trophic system.
* `RStudio_Ushio` contains the script to analyse the time series dataset of the Ushio system.


## Scripting structure

The label `RStudio` and `cmdline` in the name denote two different script structures:
* `RStudio`: corresponds to a compact implementation of NODEBNGM, meant to be executed in RStudio.
* `cmdline`: corresponds to a more modular implementation of NODEBNGM, meant to be executed via command line.

Each script folder contains the data necessary to perform the analysis, generally located in the folder `data`, as well as an `out` folder to save the output of the scripts.


