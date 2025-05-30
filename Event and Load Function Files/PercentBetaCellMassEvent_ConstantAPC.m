%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PercentBetaCellMassEvent_ConstantAPC.m
%
% Event function used in ODE solver to stop integration when beta cell mass (B)
% drops to 20% of its initial value (B_initcond), indicating critical disease threshold.
%
% Inputs:
% - ~         : time input (not used)
% - Y         : state vector where Y(5) corresponds to beta cell mass (B)
% - ~         : derivative input (not used)
% - B_initcond: initial beta cell mass
%
% Outputs:
% - value     : difference between current beta cell mass and 20% of initial mass
% - isterminal: flag to stop the solver when event occurs (1 = stop)
% - direction : zero crossing detection in both directions (0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [value, isterminal, direction] = PercentBetaCellMassEvent_ConstantAPC(~, Y, ~, B_initcond)
    B = Y(5);
    value = B - 0.2*B_initcond;  % Stop when B reaches 20% of initial value
    isterminal = 1;  % Stop integration
    direction = 0;   % Detect zero crossing in any direction
end
