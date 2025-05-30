%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% PercentBetaCellMassEvent.m
%
% Event function used with ODE solvers to stop simulation 
% when the beta cell mass (B) falls to 20% of its initial value.
%
% Inputs:
%   - ~ : time input (unused)
%   - Y : vector of system states; B is the 7th component
%   - ~ : model output or unused parameter (ignored here)
%   - B_initcond : initial beta cell mass
%
% Outputs:
%   - value : difference between current beta cell mass and 
%             20% of initial beta cell mass (zero crossing triggers event)
%   - isterminal : 1 to stop integration when event occurs
%   - direction : 0 to detect zero crossing in any direction
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [value, isterminal, direction] = PercentBetaCellMassEvent(~, Y, ~, B_initcond)
    B = Y(7);                            % Extract current beta cell mass
    value = B - 0.2*B_initcond;         % Trigger event when B reaches 20% of initial
    isterminal = 1;                     % Stop the solver when event is triggered
    direction = 0;                      % Detect zero crossing in both directions
end