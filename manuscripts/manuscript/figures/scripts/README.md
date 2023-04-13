# Scripts 

Author: Willem Bonnaffé (w.bonnaffe@gmail.com)

-----------
## Overview 

This folder contains all the scripts and data necessary to replicate the results presented in the manuscript.


-----------
## Repository structure

Each folder contains the scripts of a specific analysis of the manuscript:
* `benchmark` contains the scripts to compare the runtime and accuracy of fitting NODEs by BNGM to standard NODEs, parametric polynomial ODEs, and convergent cross mapping (CCM).
* `cmdline_HL` contains the scripts to analyse of the hare-lynx time series.
* `cmdline_TS1-3` contain the scripts to analyse the three replicated time series of the algae, flagellate, rotifer tri-trophic system.
* `RStudio_Ushio` contains the script to analyse the time series dataset of the Ushio system.


-----------
## Script structure

The label `RStudio` and `cmdline` in the name of the folders denote two different script structures:
* `RStudio`: corresponds to a compact implementation of NODEBNGM, meant to be executed in RStudio.
* `cmdline`: corresponds to a more modular implementation of NODEBNGM, meant to be executed via command line.

Each script folder contains the data necessary to perform the analysis, generally located in the folder `data`, as well as an `out` folder to save the output of the scripts.


-----------
## Installation

* `R` (v4.0.2 or later)
* Library `rEDM` (v1.14.0)


-----------
## Running the scripts

Specific instructions on how to run the scripts can be found in the `README.md` file in each folder.
