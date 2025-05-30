%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: NuRegTcellmodel_withdepletion
%
% Description:
% This function defines the system of ODEs for a regulatory T cell model 
% in the context of autoimmune diabetes, including time-dependent depletion 
% of effector and regulatory cells after a specified intervention time.
%
% Inputs:
% - t: Current time
% - Y: State vector [AL, EL, RL, AP, EP, RP, B]
% - params: Model parameters (15 values)
% - timing: Time at which depletion occurs
%
% Output:
% - DyDt: Time derivatives of the state vector
%
% Compartments:
% AL, EL, RL - Activated, Effector, Regulatory T cells in Lymph Node
% AP, EP, RP - Same cell types in Pancreas
% B         - Beta cell mass
%
% Author: [Your Name]
% Date: [Today's Date]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function DyDt = NuRegTcellmodel_withdepletion(t, Y, params, timing)

% === Unpack Model Parameters ===
phiA     = params(1);   % Migration rate of activated T cells
deltaA   = params(2);   % Death rate of activated T cells
lambdaEL = params(3);   % Activation rate of effector T cells (LN)
omegaEL  = params(4);   % Proliferation rate of EL
phiE     = params(5);   % Migration rate of EL to pancreas
deltaE   = params(6);   % Death rate of EL and EP
lambdaR  = params(7);   % Activation rate of regulatory T cells (LN)
omegaR   = params(8);   % Proliferation rate of RL
C        = params(9);   % Carrying capacity of lymph node
phiR     = params(10);  % Migration rate of RL to pancreas
deltaR   = params(11);  % Death rate of RL and RP
kappa    = params(12);  % Beta cell killing rate
r1       = params(13);  % Suppression of EL activation by RL
r2       = params(14);  % Suppression of beta cell killing by RP
nu       = params(15);  % EP activation strength

% === Relabel State Variables for Readability ===
AL = Y(1);  % Activated T cells in LN
EL = Y(2);  % Effector T cells in LN
RL = Y(3);  % Regulatory T cells in LN
AP = Y(4);  % Activated T cells in Pancreas
EP = Y(5);  % Effector T cells in Pancreas
RP = Y(6);  % Regulatory T cells in Pancreas
B  = Y(7);  % Beta cell population

% === Apply Depletion Indicator Based on Time ===
if t >= timing
    m = 0;  % Depletion has occurred
else
    m = 1;  % Pre-depletion
end

% === ODE Equations ===

% Lymph Node Dynamics
dALdt = phiA * AP - deltaA * AL;
dELdt = m * (lambdaEL * AL) / (1 + r1 * RL) ...
        + omegaEL * (1 - (EL + RL) / C) * EL ...
        - phiE * EL - deltaE * EL;
dRLdt = m * lambdaR * AL ...
        + omegaR * (1 - (EL + RL) / C) * RL ...
        - phiR * RL - deltaR * RL;

% Pancreas Dynamics
dAPdt = (nu * EP * B) / (1 + r2 * RP) - phiA * AP - deltaA * AP;
dEPdt = phiE * EL - deltaE * EP;
dRPdt = phiR * RL - deltaR * RP;
dBdt  = - (kappa * EP * B) / (1 + r2 * RP);  % Beta cell death

% === Return Derivative Vector ===
DyDt = [dALdt; dELdt; dRLdt; dAPdt; dEPdt; dRPdt; dBdt];

end
