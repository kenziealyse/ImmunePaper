# Project README



## Plotting Files

# OptimalDosingTimePlot.m
This MATLAB script simulates the effect of different immunotherapy dosing 
times on disease progression by modeling beta cell dynamics and 
plotting the relationship between effector peptide levels at dosing and 
time to critical beta cell loss.

## r1_vs_nu_heatmap.m  
Simulates and plots a heatmap of the time to 20% beta cell mass loss 
as a function of parameters \( r_1 \) and \(\nu\), visualizing 
disease progression sensitivity.

## r2_vs_nu_heatmap.m  
Simulates and plots a heatmap of time to 20% beta cell mass loss 
varying \( r_2 \) and \(\nu\), saving results for further analysis.

## TTD_vs_nu.m  
Plots the time to disease progression versus \(\nu\) for different 
initial regulatory T-cell populations, showing the impact of initial 
conditions on disease onset.

## ParameterSA.m
Performs local sensitivity analysis by varying each key model parameter 
±10% to evaluate its effect on the time to 20% beta cell mass loss, 
visualizing how parameter changes influence disease progression time.

## PhasePlanePlots.m
This script simulates and visualizes phase plane trajectories and 
beta cell mass dynamics over time for varying initial regulatory 
T-cell populations to study disease progression.

## ConstantA_vs_TTD.m
Runs simulations of the NuRegTcellmodel_constantA model
varying the parameter 'a' over a range of values to analyze its effect
on the time it takes for beta cell mass to drop to 20% of initial mass.

## plot_TTD_boundary_heatmap.m
This script loads two MATLAB data files containing simulation results
of TTD (Time to disease) as a function of parameters r1 and r2.
Part 1: Detects boundary points where TTD reaches 70 for the first 
time along each row.
Part 2: Plots a heatmap of TTD values and overlays the detected 
boundary.

## Dosing_noDepletion_diffheatmap.m
This script simulates disease progression in a T1D model under varying 
immune regulatory T-cell (RL) dosing amounts and times. For each combination 
of dose amount and timing, it computes the resulting delay (or acceleration) 
in time to 20% beta cell mass loss compared to the untreated case.

## Dosing_withAPCdepletion_diffheatmap.m
This script simulates the effects of varying both the timing and dose
of regulatory T cell therapy and APC depletion. The results are 
visualized as a heatmap showing the difference in time
to 20% beta cell mass compared to baseline disease progression.


## Function Files

## PercentBetaCellMassEvent.m 
Event function to stop ODE integration when beta cell mass falls to 
20% of its initial value, marking a critical disease threshold.

## PercentBetaCellMassEvent_ConstantAPC.m
Event function to stop ODE integration when beta cell mass falls to 
20% of its initial value, marking a critical disease threshold, for
model with constant APC.

## LoadParameters.m 
Defines and returns a vector of fixed biological and model parameters 
used in simulations.

## NuRegTcellmodel.m  
Contains the system of ODEs modeling immune cell dynamics and beta cell 
mass in the pancreas, parameterized by transition rates and regulatory 
factors.

## NuRegTcellmodel_constantA.m  
Contains the system of ODEs with constant APC modeling immune cell 
dynamics and beta cell mass in the pancreas, parameterized by 
transition rates and regulatory factors.

## NuRegTcellmodel_withdepletion.m
Contains the system of ODEs when depleting APC modeling immune cell 
dynamics and beta cell mass in the pancreas, parameterized by 
transition rates and regulatory factors.