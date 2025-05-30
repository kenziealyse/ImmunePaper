# Project README



## Plotting Files##

## r1vsNu_heatmap.m  
Simulates and plots a heatmap of the time to 20% beta cell mass loss as a function of parameters \( r_1 \) and \(\nu\), visualizing disease progression sensitivity.

## r2vsNu_heatmap.m  
Simulates and plots a heatmap of time to 20% beta cell mass loss varying \( r_2 \) and \(\nu\), saving results for further analysis.

## timeToDisease_vs_nu.m  
Plots the time to disease progression versus \(\nu\) for different initial regulatory T-cell populations, showing the impact of initial conditions on disease onset.

## ParameterSensitivityAnalysis.m
Performs local sensitivity analysis by varying each key model parameter ±10% to evaluate its effect on the time to 20% beta cell mass loss, visualizing how parameter changes influence disease progression time.



**Function Files**

**PercentBetaCellMassEvent.m**  
Event function to stop ODE integration when beta cell mass falls to 20% of its initial value, marking a critical disease threshold.

**LoadParameters.m**  
Defines and returns a vector of fixed biological and model parameters used in simulations.

**NuRegTcellmodel.m**  
Contains the system of ODEs modeling immune cell dynamics and beta cell mass in the pancreas, parameterized by transition rates and regulatory factors.