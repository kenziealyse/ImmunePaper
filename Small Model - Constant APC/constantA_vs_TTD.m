%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameter Sweep of "a" Parameter
% ----------------------------------
% This script runs simulations of the NuRegTcellmodel_constantA model
% varying the parameter 'a' over a range of values to analyze its effect
% on the time it takes for beta cell mass to drop to 20% of initial mass.
%
% Two initial conditions for regulatory T-cell population RL(0) are used:
%  - RL(0) = 10 (black dashed line)
%  - RL(0) = 1e6 (red dashed line)
%
% The results are plotted as time (years) versus 'a', and saved as a PDF.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CLEAR THE WORKSPACE
clc
clear
close all

% Load model parameters
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
nu = 0.01;

% User Defined Disease Points and Type
earlydiseasepoint = [0.2*10^(-5), 0.2*10^(-5)];
latediseasepoint = [0.8*10^(-5), 1*10^(-5)];
disease_type = 'late';

if strcmp(disease_type, 'early')
    r1 = earlydiseasepoint(1);
    r2 = earlydiseasepoint(2);
elseif strcmp(disease_type, 'late')
    r1 = latediseasepoint(1);
    r2 = latediseasepoint(2);
end

% Parameter sweep setup
N = 300;
a_vals =  linspace(0, 1e7, N); 

figurename = '../Figures\timetodiseasevsnu.pdf';

% Initial conditions for compartments
EL_initcond = 10;
RL_initcond = 10;
EP_initcond = 0;
RP_initcond = 0;
B_initcond = 1e6;

init_cond = [EL_initcond RL_initcond ...
             EP_initcond RP_initcond B_initcond]';

tspan1 = 0:0.01:70*365;

% Preallocate time vector
time = zeros(length(a_vals), 1);
  
% Simulate for initial RL = 10
for i = 1:length(a_vals)

    a = a_vals(i);

    params = [phiA deltaA lambdaEL omegaEL phiE deltaE lambdaR omegaR C ...
                  phiR deltaR kappa r1 r2 nu a]';

    % Run the Model with event detection for beta cell mass threshold
    options = odeset('Events', @(t, Y) PercentBetaCellMassEvent_ConstantAPC(t, Y, NuRegTcellmodel_constantA(t,Y, params), B_initcond));
    [T,Y] = ode23s(@(t,Y) NuRegTcellmodel_constantA(t,Y, params), tspan1, init_cond, options);

    % Save time to event (in years)
    time(i) = T(end)/365;

end

figure(1)
plot(a_vals, time, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '--')
hold on

%% Repeat simulation with increased RL(0)
RL_initcond = 1e6;

init_cond = [EL_initcond RL_initcond ...
             EP_initcond RP_initcond B_initcond]';

time = zeros(length(a_vals), 1);

for i = 1:length(a_vals)

    a = a_vals(i);

    params = [phiA deltaA lambdaEL omegaEL phiE deltaE lambdaR omegaR C ...
                  phiR deltaR kappa r1 r2 nu a]';

    options = odeset('Events', @(t, Y) PercentBetaCellMassEvent_ConstantAPC(t, Y, NuRegTcellmodel_constantA(t,Y, params), B_initcond));
    [T,Y] = ode23s(@(t,Y) NuRegTcellmodel_constantA(t,Y, params), tspan1, init_cond, options);

    time(i) = T(end)/365;

end

plot(a_vals, time, 'Color', 'r', 'LineWidth', 1.5, 'LineStyle', '--')

% Final plot formatting
ax = gca;
ax.FontSize = 21;
ylabel('Time to 20% Beta Cell Mass')
xlabel('a', 'FontSize', 30)
legend('R_L(0) = 10', 'R_L(0) = 1 \times 10^6', 'Location', 'best')
set(gcf, 'Color', 'White')
set(gcf, 'Position', [100, 300, 800, 500]);

% Save figure as PDF
set(gcf, 'Units', 'Inches');
pos = get(gcf, 'Position');
set(gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)]);
saveas(gcf, figurename);
