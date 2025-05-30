%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%
%
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function DyDt = RegTcellDosingModel_New(t, Y, params, Dose_amount, Dose_Time)

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
RLA = Y(2);
RL = Y(3);
EL = Y(4);
AP = Y(5);
EP = Y(6);
RP = Y(7);
B = Y(8);

Ueval = Dose(t, Dose_amount, Dose_Time);

if Ueval > 0
    disp(Ueval)
end

% ODE Equations
%LN
dALdt = phiA*AP - deltaA*AL;
dRLAdt = Ueval - lambdaR*AL*RLA;
dELdt = (lambdaEL*AL)/(1 + r1*RL) + omegaEL*(1 - (EL + RL)/C)*EL - phiE*EL- deltaE*EL;
dRLdt = lambdaR*AL*(1 + RLA) + (omegaR)*(1 - (EL + RL)/C)*RL - phiR*RL - deltaR*RL;

% Pancreas
dAPdt = (nu*EP*B)/(1 + r2*RP) - phiA*AP - deltaA*AP;
dEPdt = phiE*EL - deltaE*EP;
dRPdt = phiR*RL - deltaR*RP;
dBdt = -(kappa*EP*B)/(1+r2*RP);

%Output
DyDt = [dALdt; dRLAdt; dRLdt; dELdt; dAPdt; dEPdt; dRPdt; dBdt];

end