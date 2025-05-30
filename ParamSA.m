%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameter Sensitivity Analysis
%
% This script performs a local sensitivity analysis on key model parameters 
% by varying each parameter ±10% around its nominal value and simulating 
% the time to 20% beta cell mass loss. Results are plotted to visualize 
% how sensitive disease progression time is to each parameter.
%
% The model is solved using the NuRegTcellmodel ODE function, with event 
% stopping when beta cell mass reaches 20% of initial.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CLEAR THE WORKSPACE
clc
clear
close all

% Load model parameters
params = LoadParameters();

% Assign parameters for readability
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
nu = 0.01;  % Fixed value of nu

% Define disease points for regulatory parameters r1 and r2
earlydiseasepoint = [0.2*10^(-5), 0.2*10^(-5)];
latediseasepoint = [0.8*10^(-5), 1*10^(-5)];
r1 = latediseasepoint(1);
r2 = latediseasepoint(2);

% Initial conditions for compartments
AL_initcond = 0;
EL_initcond = 10;
RL_initcond = 10;
AP_initcond = 0;
EP_initcond = 0;
RP_initcond = 0;
B_initcond = 1e6;

% Time span for simulation (in days)
tspan1 = 0:0.01:70*365;

% Pack initial conditions into vector
init_cond = [AL_initcond EL_initcond RL_initcond AP_initcond ...
             EP_initcond RP_initcond B_initcond]';

% Names of parameters for labeling plots
param_names = {'\phi_A', '\delta_A', '\lambda_{E}', '\Omega_E', '\phi_E', '\delta_E',...
               '\lambda_R', '\Omega_R', 'C', '\phi_R', '\delta_R', '\kappa'};

% Original parameter vector including r1, r2, and nu
orig_params = [phiA deltaA lambdaEL omegaEL phiE deltaE lambdaR omegaR C ...
               phiR deltaR kappa r1 r2 nu]';

% Loop over each parameter to perform sensitivity analysis
for j = 1:length(param_names)
    params = orig_params;  % Reset parameters to original
    SA_upperbound = params(j) + 0.1*params(j);  % +10%
    SA_lowerbound = params(j) - 0.1*params(j);  % -10%
    % Generate 100 random samples within ±10% plus the original parameter value
    SA_param = [SA_lowerbound + (SA_upperbound - SA_lowerbound) * rand(1, 100), params(j)];
    SA_param = unique(SA_param);  % Remove duplicates
    SA_param = sort(SA_param);    % Sort ascending
    
    time = zeros(length(SA_param), 1);  % Preallocate time vector
    
    % Simulate for each parameter value
    for i = 1:length(SA_param)
        params(j) = SA_param(i);
        options = odeset('Events', @(t, Y) PercentBetaCellMassEvent(t, Y, NuRegTcellmodel(t,Y, params), B_initcond));
        [T,Y] = ode23s(@(t,Y) NuRegTcellmodel(t,Y, params), tspan1, init_cond, options);
        
        time(i) = T(end)./365;  % Convert days to years
    end
    
    % Plot results for this parameter
    subplot(3,4,j)
    plot((SA_param - orig_params(j))/orig_params(j), time, 'LineWidth', 1.3, 'Color', 'k')
    hold on
    index = find(SA_param == orig_params(j));
    orig_time = time(index);
    plot(0, orig_time, '*', 'LineWidth', 1.5, 'Color', 'r');
    xlabel(strcat(string(param_names(j)), ' fold change'))
    ylabel('Time to 20% Beta Cell Mass, years')
    ylim([50 70])
end

% Maximize figure window for better visibility
screenSize = get(0, 'ScreenSize'); % [left bottom width height]
set(gcf, 'Position', screenSize);

% Save the figure as PDF
figurename = 'Figures/paramSA.pdf';
set(gcf, 'Units', 'Inches');
pos = get(gcf, 'Position');
set(gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)]);
saveas(gcf, figurename);
