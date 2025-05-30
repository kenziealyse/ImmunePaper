# Project README



## Plotting Files

## NuRegTcellmodel.m  
Defines the full ODE system modeling immune cell dynamics and 
beta cell mass, including naive T-cell compartment transitions.


**Funtion Files**

**LoadParameters.m**  
Defines and returns the biological and model parameters used in the 
simulations.

**PercentBetaCellMassEvent.m**  
Event function that halts ODE integration when beta cell mass falls 
to 20% of its initial value, indicating critical disease progression.

**NuRegTcellmodel_constantA.m**  
Defines the system of ODEs modeling immune cell and beta cell 
dynamics with a constant source term replacing the naive T-cell 
compartment.

