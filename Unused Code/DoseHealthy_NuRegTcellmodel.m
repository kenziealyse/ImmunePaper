%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%
%
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function DyDt = DoseHealthy_NuRegTcellmodel(t, Y, params, Dose_amount, Dose_Time)

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

Ueval = Dose(t, Dose_amount, Dose_Time);

if Ueval > 0
    disp(Ueval)
end

% ODE Equations
%LN
dALdt = phiA*AP - deltaA*AL;
dELdt = (lambdaEL*AL)/(1 + r1*RL) + omegaEL*(1 - (EL + RL)/C)*EL - phiE*EL- deltaE*EL;
dRLdt = lambdaR*AL + (omegaR)*(1 - (EL + RL)/C)*RL - phiR*RL - deltaR*RL + Ueval;

% Pancreas
dAPdt = (nu*EP*B)/(1 + r2*RP) - phiA*AP - deltaA*AP;
dEPdt = phiE*EL - deltaE*EP;
dRPdt = phiR*RL - deltaR*RP;
dBdt = -(kappa*EP*B)/(1+r2*RP);

%Output
DyDt = [dALdt; dELdt; dRLdt; dAPdt; dEPdt; dRPdt; dBdt];

end