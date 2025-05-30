%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function defines the system of ODEs modeling the dynamics of 
% regulatory and effector T-cells, antigen-presenting cells, and beta cell mass.
%
% Inputs:
%   - ~      : time variable (not used explicitly)
%   - Y       : state vector [AL, EL, RL, AP, EP, RP, B]
%   - params  : vector of model parameters [phiA, deltaA, lambdaEL, omegaEL, phiE,
%               deltaE, lambdaR, omegaR, C, phiR, deltaR, kappa, r1, r2, nu]
%
% Outputs:
%   - DyDt    : derivatives of the state variables (rate of change)
%
% The model includes dynamics in the lymph nodes (AL, EL, RL), pancreas (AP, EP, RP),
% and beta cell population (B), accounting for proliferation, activation, decay, 
% and suppression mechanisms modulated by parameters r1 and r2.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function DyDt = NuRegTcellmodel(~, Y, params)

% Paramter Values
phiA = params(1);
deltaA = params(2);
lambdaEL = params(3);
omegaEL = params(4);
phiE = params(5);
deltaE = params(6);
lambdaR = params(7);
omegaR = params(8);
C = params(9);
phiR = params(10);
deltaR = params(11);
kappa = params(12);
r1 = params(13);
r2 = params(14);
nu = params(15);

% Relabel Compartments to Easily Keep Track
AL = Y(1);
EL = Y(2);
RL = Y(3);
AP = Y(4);
EP = Y(5);
RP = Y(6);
B = Y(7);

% ODE Equations
%LN
dALdt = phiA*AP - deltaA*AL;
dELdt = (lambdaEL*AL)/(1 + r1*RL) + omegaEL*(1 - (EL + RL)/C)*EL - phiE*EL- deltaE*EL;
dRLdt = lambdaR*AL + (omegaR)*(1 - (EL + RL)/C)*RL - phiR*RL - deltaR*RL;

% Pancreas
dAPdt = (nu*EP*B)/(1 + r2*RP) - phiA*AP - deltaA*AP;
dEPdt = phiE*EL - deltaE*EP;
dRPdt = phiR*RL - deltaR*RP;
dBdt = -(kappa*EP*B)/(1+r2*RP);


%Output
DyDt = [dALdt; dELdt; dRLdt; dAPdt; dEPdt; dRPdt; dBdt];

end