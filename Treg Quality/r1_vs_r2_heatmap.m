%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script Name: r1_vs_r2_heatmap.m
%
% Description:
% This script evaluates the effect of varying immune regulation parameters 
% (r₁ and r₂) on disease progression in a type 1 diabetes (T1D) model.
% It simulates the time to 20% beta cell mass for a grid of r₁ and r₂ values 
% and visualizes the results in a heatmap.
%
% Output:
% - Heatmap figure showing time to disease (in years) as a function of r₁ and r₂
% - Saved .mat file containing simulation data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CLEAR THE WORKSPACE
clc
clear
close all

% Time span: simulate for up to 70 years
tspan = 0:1:70*365;

% Load base model parameters
params = LoadParameters();

% Unpack fixed model parameters
phiA = params(1); deltaA = params(2); lambdaEL = params(3);
omegaEL = params(4); phiE = params(5); deltaE = params(6);
lambdaR = params(7); omegaR = params(8); C = params(9);
phiR = params(10); deltaR = params(11); kappa = params(15);
r1 = params(16);  % Overwritten in loop
nu = 0.01;        % Regulatory strength

% Initial Conditions
AL_initcond = 0;
EL_initcond = 10;
RL_initcond = 1e6;    % High RL initial condition
AP_initcond = 0;
EP_initcond = 0;
RP_initcond = 0;
B_initcond = 1e6;

init_cond = [AL_initcond EL_initcond RL_initcond AP_initcond ...
             EP_initcond RP_initcond B_initcond]';

% Output figure path
figurename = ['..\Figures\highRL0r1vsr2heatmapnuis', num2str(nu), '.pdf'];

% Define r₁ and r₂ parameter space
N = 300;
r1_vals = linspace(0, 2e-4, N);
r2_vals = linspace(0, 2e-4, N);

% Preallocate matrix to store time to disease
T_vals = zeros(length(r2_vals), length(r1_vals));

% Loop over r₁ and r₂ combinations
for i = 1:length(r1_vals)
    r1 = r1_vals(i);

    for j = 1:length(r2_vals)
        r2 = r2_vals(j);

        % Update parameter vector with current r₁ and r₂
        params = [phiA deltaA lambdaEL omegaEL phiE deltaE lambdaR omegaR C ...
                  phiR deltaR kappa r1 r2 nu]';

        % Solve ODE system with event to stop at 20% beta cell mass
        options = odeset('Events', @(t, Y) ...
            PercentBetaCellMassEvent(t, Y, NuRegTcellmodel(t,Y, params), B_initcond));
        [T, Y] = ode23s(@(t, Y) NuRegTcellmodel(t, Y, params), tspan, init_cond, options);

        % Record time to event (converted to years)
        T_vals(j, i) = T(end) / 365;
    end
end

% Save results to .mat file
results.r2vals = r2_vals;
results.r1_vals = r1_vals;
results.TTD = T_vals;
save(['Mat Files/highRL0r1vsr2heatmap', num2str(nu), '_vals.mat'], 'results')

% -----------------------------
% Plotting Heatmap
% -----------------------------
figure(1)
imagesc(r1_vals, r2_vals, T_vals);
set(gca, 'YDir', 'normal')  % Standard orientation
xlim([0 3e-5])
ylim([0 3e-5])
xlabel('r_1', 'FontSize', 30);
ylabel('r_2', 'FontSize', 30);
title('Time to 20% Beta Cell Mass', 'FontSize', 21);
set(gca, 'YScale', 'log')   % Log scale for r₂ axis
ax = gca;
ax.FontSize = 21;

% Define custom RGB colormap
numColors = 100;
redValues   = [linspace(1, 0, numColors/2), zeros(1, numColors/2)];
greenValues = [zeros(1, numColors/2), linspace(0, 1, numColors/2)];
blueValues  = [linspace(0, 1, numColors/2), linspace(1, 0, numColors/2)];
customColorMap = [redValues', greenValues', blueValues']; 
colormap(gca, customColorMap);

% Colorbar and limits
colorbar;
clim([0, 70]);  % Maximum display range in years

% Set figure size and save
set(gcf, 'Position', [100, 300, 800, 500]);
set(gcf, 'Color', 'White');
set(gcf, 'Units', 'Inches');
pos = get(gcf, 'Position');
set(gcf, 'PaperPositionMode', 'Auto', ...
         'PaperUnits', 'Inches', ...
         'PaperSize', [pos(3), pos(4)]);
saveas(gcf, figurename);