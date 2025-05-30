%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%
%
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function DyDt = SmallModel2_ND(~, Y, params)

% Paramter Values
p1 = params(1);
p2 = params(2);
p3 = params(3);


% Relabel Compartments to Easily Keep Track
EL = Y(1);
RL = Y(2);
AP = Y(3);
B = Y(4);

% ODE Equations
%LN
dELdt = AP/(1 + p3*RL) + (1 - (EL + RL))*EL;
dRLdt = AP + (1 - (EL + RL))*RL;

% Pancreas
dAPdt = EL*B - (p2)*AP;
dBdt = -(p1)*EL*B;

%Output
DyDt = [dELdt; dRLdt; dAPdt; dBdt];

end