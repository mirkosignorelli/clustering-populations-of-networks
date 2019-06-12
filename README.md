# clustering-populations-of-networks

Data and code to reproduce the simulations in:

Signorelli, M., Wit, E. C., Model-based clustering for populations of networks. arXiv preprint: arXiv:1806.00225.

The code to run each simulation is divided into scripts sequentially numbered:
- script 1 generates the data
- script 2 runs the Expectation-Maximization algorithm for all subcases and repetitions considered (this may in some cases be time consuming, so we provide an implementation that allows to easily perform parallel computing)
- script 3 further processes the output of script 2 to generate the results presented in the paper (with the exception of simulation J, where first model selection is performed in script 3 and then the results are produced in script 4).

Note that scripts 1A, 1B, ..., 2A, 2B, ... require to source some scripts whose name start with 0-; these scripts can be found in the folder "functions to source".
