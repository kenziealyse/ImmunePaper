%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script simulates a regulatory T-cell model over a grid of parameter values
% for r₁ (Reg T cell suppression rate) and ν (activation rate of effector cells by APCs).
%
% For each (r₁, ν) pair, it computes the time until beta cell mass declines to 20%
% of its initial value and visualizes the result in a heatmap.
%
% Model parameters are loaded using LoadParameters.m, and ODEs are solved using ode23s.
% A custom colormap enhances interpretability of the heatmap.
%
% Output:
% - A PDF heatmap showing time to 20% beta cell mass across the (r₁, ν) parameter space.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Clear workspace
clc
clear
close all

%% Set time span and model parameters
tspan = 0:.1:70*365;  % simulate for 70 years (in days)
params = LoadParameters();

% Extract specific parameters from vector
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
r2 = 0;  % keep r2 fixed at 0

%% Initial conditions
AL_initcond = 0;
EL_initcond = 10;
RL_initcond = 10;
AP_initcond = 0;
EP_initcond = 0;
RP_initcond = 0;
B_initcond = 1e6;

init_cond = [AL_initcond EL_initcond RL_initcond AP_initcond ...
             EP_initcond RP_initcond B_initcond]';

%% Define r₁ and ν vectors to sweep
N = 300;
r1_vals = linspace(0, 0.6e-3, N);
nu_vals = linspace(0.25e-5, 2500e-5, N); 

T_vals = zeros(length(nu_vals), length(r1_vals));  % preallocate

%% Simulate model over parameter grid
for i = 1:length(r1_vals)
    r1 = r1_vals(i);
    for j = 1:length(nu_vals)
        nu = nu_vals(j);

        % Update parameter vector for this (r1, nu) pair
        params = [phiA deltaA lambdaEL omegaEL phiE deltaE lambdaR omegaR C ...
                  phiR deltaR kappa r1 r2 nu]';

        % Define event to stop when beta cell mass drops below 20%
        options = odeset('Events', @(t, Y) PercentBetaCellMassEvent(t, Y, NuRegTcellmodel(t,Y, params), B_initcond));

        % Run ODE model
        [T, Y] = ode23s(@(t,Y) NuRegTcellmodel(t,Y, params), tspan, init_cond, options);

        % Save final time (in years)
        T_vals(j,i) = T(end)/365;
    end  
end

%% Save results in .mat
results.r1vals = r1_vals;
results.nuvals = nu_vals;
results.TTD = T_vals;
save('Mat Files/NuHeatmapr1_vals.mat', 'results')

%% Create heatmap plot
figure(1)
imagesc(r1_vals, nu_vals, T_vals);
hold on
xlim([0 max(r1_vals)])
ylim([0 0.5e-3])
xlabel('r_1', 'FontSize', 17)
ylabel('\nu', 'FontSize', 17)
set(gca, 'YDir', 'normal')
title('Time to 20% Beta Cell Mass')
ax = gca;
ax.FontSize = 21;

%% Custom colormap setup
numColors = 100;
redValues = [linspace(1, 0, numColors/2), zeros(1, numColors/2)];
greenValues = [zeros(1, numColors/2), linspace(0, 1, numColors/2)];
blueValues = [linspace(0, 1, numColors/2), linspace(1, 0, numColors/2)];
customColorMap = [redValues', greenValues', blueValues']; 
colormap(gca, customColorMap);
colorbar;
clim([0, 70]);  % color scale in years

% Set figure size
set(gcf, 'Position', [100, 300, 800, 500]);

%% Save figure to file
figurename = '../Figures/r1vsNuheatmapDynamics.pdf';
set(gcf, 'Units', 'Inches');
pos = get(gcf, 'Position');
set(gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)]);
saveas(gcf, figurename);
