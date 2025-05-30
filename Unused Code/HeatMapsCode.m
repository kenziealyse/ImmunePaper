% Clear the workspace
clc
clear
close all

heatmap = 'r1vsnu';

switch heatmap

    case 'r1vsnu'

    


    case 'r2vsnu'

    case 'r1vsr2lownu'

    case 'r1vsr2highnu'

end


function [value, isterminal, direction] = events(~, Y, ~, B_initcond)
    B = Y(7);
    value = B - 0.2*B_initcond;  % Condition: stop when B is close to 0.2*initcondB
    isterminal = 1;  % Stop the integration
    direction = -1;   % Detect all zeros (both rising and falling)
end