%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%
%
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function DyDt = ToyModel2_ND(t, Y, params)

% Paramter Values
deltaA = params(2);
omegaR = params(8);
C = params(9);
kappa = params(12);
r1 = params(13);
r2 = params(14);

% Relabel Compartments to Easily Keep Track
EL = Y(1);
RL = Y(2);
AP = Y(3);
B = Y(4);

if t/(omegaR*365) > 5
    m = 0;
else
    m = 1;
end
m=1;
% ODE Equations
%LN
dELdt = m*(AP)/(1 + r1*C*RL) + (1 - (EL + RL))*EL;
dRLdt = m*AP + (1 - (EL + RL))*RL;

% Pancreas
dAPdt = (EL*B)/(1 + r2*C*RL) - deltaA/omegaR*AP;
dBdt = -(kappa*C*EL*B)/(omegaR*(1+r2*C*RL));

%Output
DyDt = [dELdt; dRLdt; dAPdt; dBdt];

end