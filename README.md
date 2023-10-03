# Replication files for "Combining forecasts with mean absolute error scoring rule"

## Content
This repository contains codes, data, and other supplementary materials by [Felix Chan](mailto:F.Chan@curtin.edu.au) and [Laurent Pauwels](mailto:laurent.pauwels@sydney.edu.au) for the paper entitled 

	Chan, F. and Pauwels, L. "Optimal forecast combination with mean absolute error loss", 2023. 

If you are using the code, please cite the paper.

## Software requirements
The software used for replications are:

 - MATLAB Version: 9.14.0.2206163 (R2023a)
 - Julia 1.6.1
 
For the full Julia requirements check *Requirements.txt* and each Jupyter Notebook for the specific packages used. The simulations and empirical illustrations presented in the paper are conducted with MATLAB and Julia. The computer code is available in both languages. The Julia Language code is presented in a Jupyter Notebook. There may be variations in the results because of differences in the random number generator and the solver (optimizer) in the two software.
	
## code
  
### julia
The Julia Language code is presented in a Jupyter Notebook.

1. `TheoryGraphs.ipynb`: Jupyter Notebook to construct the graphs featured in Figure 1 of Section 2 of the paper. The framework with which those graphs are constructed is explained in the Jupyter Notebook. 

2. `SimulationsJulia.ipynb`: Jupyter Notebook with Julia Language to conduct the simulations found in Section 4. 


3. `EmpiricsJulia.ipynb`: Jupyter Notebook with Julia Language to conduct the empirical study found in Section 5. 

### matlab

1. `Simulations.m`: MATLAB script that conducts the simulations found in Section 4 of the paper. 

2. `Empirics.m`: MATLAB script that conducts the empirical study found in Section 5 of the paper. 

3. `maeloss.m`: MATLAB function that computes mean absolute error loss. It is used with `Simulations.m` and `Empirics.m`.

## data

### empirical
The datasets for the empirical illustrations in `.csv` and MATLAB `.mat` formats are:
 
- `InflationData`: the inflation rate data.
- `GrowthData`: the real growth rate data.
- `UnempData`: the unemployment rate data.

We focus on the expert-predicted rate of inflation, real GDP growth, and unemployment for the European Union for the coming year (one year ahead). The data are from the following sources:

- Source: ECB Survey of Professional Forecasters is available at <https://www.ecb.europa.eu/stats/ecb_surveys/survey_of_professional_forecasters/html/index.en.html>.

### simulated

The simulated data which replicates the simulation results are:

- `sim_t3Set1.csv` and `sim_snSet1.csv` or `SN_Set1.mat`.
- `sim_t3Set2.csv` and `sim_snSet2.csv` or `SN_Set2.mat` .

These files contain the exact values for `b` (Skewness Parameter) and `Sigma` (Random covariance matrix) for 2 different sets of optimal weights, `a`.


## paper
    This directory contains the supplementary materials for the paper. 

    
