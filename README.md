# Model-based clustering for populations of networks

This folder contains R scripts related to the publication of Signorelli, M., Wit, E. C. (2020). Model-based clustering for populations of networks. *Statistical Modelling*, 20 (1).
You can read the paper (with open access) here: https://journals.sagepub.com/doi/full/10.1177/1471082X19871128

## About this repository
This repository contains the data and code to reproduce the simulations presented in Signorelli and Wit (2020).

The code to run each simulation is divided into scripts sequentially numbered:
- script 1 generates the data
- script 2 runs the Expectation-Maximization algorithm for all subcases and repetitions considered (this may in some cases be time consuming, so we provide an implementation that allows to easily perform parallel computing)
- script 3 further processes the output of script 2 to generate the results presented in the paper (with the exception of simulation J, where first model selection is performed in script 3 and then the results are produced in script 4).

Note that scripts 1A, 1B, ..., 2A, 2B, ... require to source some scripts whose name start with 0-; these scripts can be found in the folder "functions to source".
