%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script simulates a regulatory T-cell model across a grid of parameter values
% for r₂ (alternative Reg T cell suppression mechanism) and ν (activation rate of effector cells by APCs).
%
% For each (r₂, ν) pair, it calculates the time until beta cell mass declines to 20%
% of its initial value. The results are saved in a .mat file and visualized as a heatmap.
%
% Output:
% - Heatmap figure showing time to 20% beta cell mass across (r₂, ν) space
% - Results saved to 'Mat Files/NuHeatmapr2_vals.mat'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Clear workspace and setup
clc
clear
close all

% Simulation time span (in days)
tspan = 0:1:70*365;

% Load base parameters
params = LoadParameters();

% Extract specific parameters
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
r1 = 0;  % keep r1 fixed at 0

% Initial conditions
AL_initcond = 0;
EL_initcond = 10;
RL_initcond = 10;
AP_initcond = 0;
EP_initcond = 0;
RP_initcond = 0;
B_initcond = 1e6;

init_cond = [AL_initcond EL_initcond RL_initcond AP_initcond ...
             EP_initcond RP_initcond B_initcond]';

% Output filename
figurename = '../Figures/r2vsNuheatmapDynamics.pdf';

%% Parameter sweep setup
N = 300;
r2_vals = linspace(0, 3e-3, N);
nu_vals = linspace(0.25e-5, 2500e-5, N);

% Preallocate result matrix
T_vals = zeros(length(nu_vals), length(r2_vals));

%% Run simulation grid
for i = 1:length(r2_vals)
    r2 = r2_vals(i)

    for j = 1:length(nu_vals)
        nu = nu_vals(j);

        % Update parameter vector
        params = [phiA deltaA lambdaEL omegaEL phiE deltaE lambdaR omegaR C ...
                  phiR deltaR kappa r1 r2 nu]';

        % Solve ODEs with event to stop at 20% beta cell mass
        options = odeset('Events', @(t, Y) PercentBetaCellMassEvent(t, Y, NuRegTcellmodel(t,Y, params), B_initcond));
        [T, Y] = ode23s(@(t,Y) NuRegTcellmodel(t,Y, params), tspan, init_cond, options);

        % Store time to depletion (in years)
        T_vals(j, i) = T(end) / 365;
    end
end

%% Save results
results.r2vals = r2_vals;
results.nuvals = nu_vals;
results.TTD = T_vals;
save('Mat Files/NuHeatmapr2_vals.mat', 'results');

%% Plot heatmap
figure(1)
imagesc(r2_vals, nu_vals, T_vals);
xlim([0 max(r2_vals)])
ylim([0 0.5e-3])
ax = gca;
ax.FontSize = 21;
xlabel('r_2', 'FontSize', 30);
ylabel('\nu', 'FontSize', 30);
set(gca, 'YDir', 'normal');
title('Time to 20% Beta Cell Mass', 'FontSize', 21);

% Define custom colormap
numColors = 100;
redValues = [linspace(1, 0, numColors/2), zeros(1, numColors/2)];
greenValues = [zeros(1, numColors/2), linspace(0, 1, numColors/2)];
blueValues = [linspace(0, 1, numColors/2), linspace(1, 0, numColors/2)];
customColorMap = [redValues', greenValues', blueValues'];

% Apply colormap
colormap(gca, customColorMap);
colorbar;
clim([0, 70]);  % color scale in years

set(gcf, 'Position', [100, 300, 800, 500]);

%% Save figure to file
set(gcf, 'Units', 'Inches');
pos = get(gcf, 'Position');
set(gcf, 'PaperPositionMode', 'Auto', ...
         'PaperUnits', 'Inches', ...
         'PaperSize', [pos(3), pos(4)]);
saveas(gcf, figurename);
