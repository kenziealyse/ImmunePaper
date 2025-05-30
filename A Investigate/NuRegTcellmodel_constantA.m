%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%
%
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function DyDt = NuRegTcellmodel_constantA(~, Y, params)

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
a = params(16);

% Relabel Compartments to Easily Keep Track
EL = Y(1);
RL = Y(2);
EP = Y(3);
RP = Y(4);
B = Y(5);

% ODE Equations
%LN
dELdt = (lambdaEL*a)/(1 + r1*RL) + omegaEL*(1 - (EL + RL)/C)*EL - phiE*EL- deltaE*EL;
dRLdt = lambdaR*a + (omegaR)*(1 - (EL + RL)/C)*RL - phiR*RL - deltaR*RL;

% Pancreas
dEPdt = phiE*EL - deltaE*EP;
dRPdt = phiR*RL - deltaR*RP;
dBdt = -(kappa*EP*B)/(1+r2*RP);


%Output
DyDt = [dELdt; dRLdt; dEPdt; dRPdt; dBdt];

end