%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%
%
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function DyDt = RegTcellModel(~, Y, params)

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
sigma = params(12);
alpha = params(13);
deltaP = params(14);
kappa = params(15);
r1 = params(16);
r2 = params(17);

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
dAPdt = (sigma*alpha*EP*B)/(deltaP*(1 + r2*RP)) - phiA*AP - deltaA*AP;
dEPdt = phiE*EL - deltaE*EP;
dRPdt = phiR*RL - deltaR*RP;
dBdt = -(kappa*EP*B)/(1+r2*RP);

%Output
DyDt = [dALdt; dELdt; dRLdt; dAPdt; dEPdt; dRPdt; dBdt];

end