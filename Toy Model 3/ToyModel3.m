%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%
%
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function DyDt = SmallModel1(~, Y, params)

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
EL = Y(1);
RL = Y(2);
AP = Y(3);
B = Y(4);

% ODE Equations
%LN
dELdt = (lambdaEL*AP)/(1 + r1*RL) + omegaEL*(1 - (EL + RL)/C)*EL;
dRLdt = lambdaR*AP + (omegaR)*(1 - (EL + RL)/C)*RL;

% Pancreas
dAPdt = (nu*EL*B) - deltaA*AP;
dBdt = -(kappa*EL*B);

%Output
DyDt = [dELdt; dRLdt; dAPdt; dBdt];

end