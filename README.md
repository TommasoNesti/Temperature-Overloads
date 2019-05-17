# Computer code for the paper "Temperature Overloads in Power Grids Under Uncertainty: a Large Deviations Approach"

## Summary  
This repository contains the code for reproducing the Figures in the paper  

[1] "Temperature Overloads in Power Grids Under Uncertainty: a Large Deviations Approach", Tommaso Nesti, Jayakrishnan Nair 
     and Bert Zwart, accepted for publication in IEEE Transactions on Control of Network Systems.

## Prerequisites  
* The script for Figure 1 requires Mathematica (https://www.wolfram.com/mathematica/)
* The script for Figure 2 requires MATLAB (https://www.mathworks.com/products/matlab.html) and MATPOWER (http://www.pserc.cornell.edu/matpower/) 

### Instructions
No installation is required. Simply download and run the scripts.
* 'generate_fig1.nb' is a Mathematica notebook used to generate Fig 1 in [1].
* 'generate_fig2.m' is a Matlab script used to generate Fig 2 in [1]. Change line 23 to cd your MATPOWER folder.
* 'solve_uncertainty_aware_OPF.m' is a Malatb function called within the script 'generate_fig2.m', and relies on the MATPOWER functions 'makePTDF' and 'rundcopf'.


