%% Event
function [value, isterminal, direction] = PercentBetaCellMassEvent(~, Y, ~, B_initcond)
    B = Y(5);
    value = B - 0.2*B_initcond;  % Condition: stop when dydt is close to zero
    isterminal = 1;  % Stop the integration
    direction = 0;   % Detect all zeros (both rising and falling)
end