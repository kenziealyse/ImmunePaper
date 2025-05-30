% CLEAR THE WORKSPACE
clc
clear
close all

params = LoadParameters();

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
kappa = params(15);

earlydiseasepoint = [0.2*10^(-5), 0.2*10^(-5)];
latediseasepoint = [0.8*10^(-5), 1*10^(-5)];

r1 = latediseasepoint(1);
r2 = latediseasepoint(2);

N = 5;

% Define the range of initial conditions for EL, RL, and nu
EL_range = [1e1, 1e2, 1e3, 1e4, 1e5, 1e6];%linspace(1, 10, N);  % Example range for EL
RL_range = [1e1, 1e2, 1e3, 1e4, 1e5, 1e6];%linspace(1, 1e6, N);  % Example range for RL
nu_vals =  linspace(1*10^(-5), 2500*10^(-5), N); 

% Initialize the result matrix to store the ratios (EL/RL) for each nu
ratio_vals = zeros(length(EL_range), length(RL_range), length(nu_vals));  % 3D matrix for EL/RL ratios

AP_initcond = 0;
B_initcond = 1*10^6;

tspan = 0:.1:70*365;

% Create a meshgrid for EL and RL to represent all combinations
[EL_matrix, RL_matrix] = meshgrid(EL_range, RL_range);

% Calculate the EL/RL ratio for each combination
for i = 1:numel(EL_matrix)
    EL_RL_ratios(i) = EL_matrix(i) / RL_matrix(i);  % Store the EL/RL ratio
end

for j = 1:length(nu_vals)  
    nu = nu_vals(j);

    for i = 1:numel(EL_matrix)
        EL_initcond = EL_matrix(i);
        RL_initcond = RL_matrix(i);

        init_cond = [EL_initcond RL_initcond AP_initcond B_initcond 0]';

        params = [phiA deltaA lambdaEL omegaEL phiE deltaE lambdaR omegaR C ...
                  phiR deltaR kappa r1 r2 nu]';

        % Run the Model
        options = odeset('Events', @(t, Y) RatioEvent(t, Y, kappa));
        [T,Y] = ode23s(@(t,Y) SmallModel2(t,Y, params), tspan, init_cond, options);

        % Store the ratio in the ratio_vals matrix (3D matrix: [EL, RL, nu])
        T_vals(i, j) = T(end)./365;  

    end
end

% Plot the heatmap
figure;
imagesc(nu_vals, EL_RL_ratios, T_vals);  % Create heatmap
colorbar;  % Display color scale
xlabel('nu');  % Label for x-axis (nu)
ylabel('EL/RL');  % Label for y-axis (EL/RL ratio)
title('Heatmap of Model Results for Varying EL/RL Ratio and nu');
set(gca, 'YDir', 'normal');  % Flip y-axis for correct visualization

% Define the number of colors
numColors = 100;

% Define transitions for red, green, and blue
redValues = [linspace(1, 0, numColors/2), zeros(1, numColors/2)]; % Red fades to zero
greenValues = [zeros(1, numColors/2), linspace(0, 1, numColors/2)]; % Green grows
blueValues = [linspace(0, 1, numColors/2), linspace(1, 0, numColors/2)]; % Blue in middle
ax = gca;
ax.FontSize = 21;
% Combine into a custom colormap
customColorMap = [redValues', greenValues', blueValues']; 

% Apply the colormap
colormap(gca, customColorMap);

% Show color bar with custom tick labels
colorbar;
clim([0, 70]);  % Set color axis limits based on data range
hold on

set(gcf, 'Position', [100, 300, 800, 500]);


function [value, isterminal, direction] = RatioEvent(~,Y, kappa)
    ratio = Y(5);
    value = ratio - log(5)/kappa;  % Condition: stop when dydt is close to zero
    % B = Y(4);
    % value = B - 0.2.*1*10^6;
    isterminal = 1;  % Stop the integration
    direction = 0;   % Detect all zeros (both rising and falling)
end