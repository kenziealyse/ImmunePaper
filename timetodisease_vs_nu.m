%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% timetodisease_vs_nu.m
%
% This script simulates the time for beta cell mass to reduce 
% to 20% under varying values of parameter \nu, representing 
% a biological rate in the T-cell regulatory model.
%
% Two scenarios are compared:
%   1) Initial regulatory T-cell count RL(0) = 10
%   2) Increased initial regulatory T-cell count RL(0) = 1e6
%
% The model integrates ODEs defined in NuRegTcellmodel.m and 
% uses event detection to stop simulation when beta cell mass 
% reaches 20% of its initial value.
%
% Outputs a plot of time-to-disease vs. \nu for both scenarios 
% and saves the figure as a PDF.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CLEAR THE WORKSPACE
clc
clear
close all

% Load model parameters
params = LoadParameters();

% Extract individual parameters for clarity
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

% Define disease points for parameters r1 and r2 (e.g. late disease state)
earlydiseasepoint = [0.2*10^(-5), 0.2*10^(-5)];
latediseasepoint = [0.8*10^(-5), 1*10^(-5)];

r1 = latediseasepoint(1);
r2 = latediseasepoint(2);

% Set number of samples and range for \nu parameter
N = 100;
nu_vals = linspace(1*10^(-5), 2500*10^(-5), N);

% Figure save location
figurename = 'Figures\timetodiseasevsnu.pdf';

% Initial conditions for compartments: [AL EL RL AP EP RP B]
AL_initcond = 0;
EL_initcond = 10;
RL_initcond = 10;       % Scenario 1 initial RL(0)
AP_initcond = 0;
EP_initcond = 0;
RP_initcond = 0;
B_initcond = 1e6;

init_cond = [AL_initcond EL_initcond RL_initcond AP_initcond EP_initcond RP_initcond B_initcond]';

tspan1 = 0:0.01:70*365; % Simulation time span (days to years)

% Preallocate array for storing time to 20% beta cell mass
time = zeros(length(nu_vals), 1);

% Loop over all nu values for scenario 1 (RL(0) = 10)
for i = 1:length(nu_vals)
    nu = nu_vals(i);
    params = [phiA deltaA lambdaEL omegaEL phiE deltaE lambdaR omegaR C ...
              phiR deltaR kappa r1 r2 nu]';

    % Solve ODE with event detection for beta cell mass threshold
    options = odeset('Events', @(t, Y) PercentBetaCellMassEvent(t, Y, NuRegTcellmodel(t,Y, params), B_initcond));
    [T,Y] = ode23s(@(t,Y) NuRegTcellmodel(t,Y, params), tspan1, init_cond, options);

    % Record time to event (convert from days to years)
    time(i) = T(end) / 365;
end

% Plot results for scenario 1
figure(1)
plot(nu_vals, time, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '--')
hold on

% Scenario 2: Increase initial regulatory T-cell count RL(0)
RL_initcond = 1e6;
init_cond = [AL_initcond EL_initcond RL_initcond AP_initcond EP_initcond RP_initcond B_initcond]';

% Preallocate again for scenario 2
time = zeros(length(nu_vals), 1);

% Loop over nu values for scenario 2 (RL(0) = 1e6)
for i = 1:length(nu_vals)
    nu = nu_vals(i);
    params = [phiA deltaA lambdaEL omegaEL phiE deltaE lambdaR omegaR C ...
              phiR deltaR kappa r1 r2 nu]';

    options = odeset('Events', @(t, Y) PercentBetaCellMassEvent(t, Y, NuRegTcellmodel(t,Y, params), B_initcond));
    [T,Y] = ode23s(@(t,Y) NuRegTcellmodel(t,Y, params), tspan1, init_cond, options);

    time(i) = T(end) / 365;
end

% Plot results for scenario 2
plot(nu_vals, time, 'Color', 'r', 'LineWidth', 1.5, 'LineStyle', '--')

% Formatting the plot
ax = gca;
ax.FontSize = 21;
ylabel('Time to 20% Beta Cell Mass (years)', 'FontSize', 20)
xlabel('\nu', 'FontSize', 30)
legend('R_L(0) = 10', 'R_L(0) = 1 \times 10^6', 'Location', 'best')

set(gcf, 'Position', [100, 300, 800, 500]);

% Save the figure as a PDF
set(gcf, 'Units', 'Inches');
pos = get(gcf, 'Position');
set(gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)]);
saveas(gcf, figurename);