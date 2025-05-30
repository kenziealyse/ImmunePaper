 
% CLEAR THE WORKSPACE
clc
clear
close all

% Set Model Conditions
tspan = 0:.1:150*365;

% Load Parameter Values
params = LoadParameters();

% Specify Parameter Values
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

% Specify Initial Conditions
AL_initcond = 0;
EL_initcond = 10;
RL_initcond = 10; %1e6;
AP_initcond = 0;
EP_initcond = 0;
RP_initcond = 0;
B_initcond = 1*10^6;

% Save Initital Conditions in a Vector
init_cond = [AL_initcond EL_initcond RL_initcond AP_initcond ...
             EP_initcond RP_initcond B_initcond]';

% Set Scaling Factor of EL and RL for Plotting
scale = 1e7;

% Define kE and kR
kR = 1/((omegaR)/(C*(omegaR - (phiR+deltaR))));
kE = 1/((omegaEL)/(C*(omegaEL - (phiE+deltaE))));

% Define r1 and r2 Values for  early and late disease and high and low nu
r1_vals = [0.2*10^(-5), 0.8*10^(-5)];
r2_vals = [0.2*10^(-5), 1*10^(-5)];
nu_vals = [0.01, 1*10^(-5)];

% Specify Figure Name
figurename = ['Figures\ELvsRLandALplot_latedisease_RLis', num2str(RL_initcond), '.pdf'];

%% Late Disease, High and Low Nu

% Specify Parameter Values
r1 = r1_vals(2);
r2 = r2_vals(2);
nu = nu_vals(1);

% Save parameters in a vector
params = [phiA deltaA lambdaEL omegaEL phiE deltaE lambdaR omegaR C ...
      phiR deltaR kappa r1 r2 nu]';

% Run the Model
options = odeset('Events', @(t, Y) events(t, Y, NuRegTcellmodel(t,Y, params), B_initcond));
[T,Y] = ode23s(@(t,Y) NuRegTcellmodel(t,Y, params), tspan, init_cond, options);

% Relabel to easily keep track of desired compartments
AL = Y(:,1);
EL = Y(:,2);
RL = Y(:,3);

disp(T(end)./365)

% Plot EL vs RL for high nu
figure(1)
subplot(1,2,1)
plot(EL./scale, RL./scale, 'LineWidth', 2, 'Color', 'b', 'LineStyle', '-')
xlabel('E_L \times 10^7', 'FontSize', 17)
ylabel('R_L \times 10^7', 'FontSize', 17)
% ylim([1.5 3])
hold on

figure(2)
subplot(2,2,1)
rate1 = lambdaR.*AL;
rate2 = (omegaR).*(1 - (EL + RL)/C).*RL - phiR.*RL - deltaR.*RL;
plot(T./365, rate1./scale, 'LineWidth', 2, 'Color', 'b', 'LineStyle', '--')
hold on
plot(T./365, rate2./scale, 'LineWidth', 2, 'Color', 'r', 'LineStyle', '-')
% legend('First Rate for RL - with AL', 'Second Rate for RL - with carying capicity')
set(gca, 'Xscale', 'log')
title('low RL, high nu')
ylim([-5 10])


% Plotting Low nu Value
nu = nu_vals(2);

% Save parameters in a vector
params = [phiA deltaA lambdaEL omegaEL phiE deltaE lambdaR omegaR C ...
      phiR deltaR kappa r1 r2 nu]';

% Run the Model
options = odeset('Events', @(t, Y) events(t, Y, NuRegTcellmodel(t,Y, params), B_initcond));
[T,Y] = ode23s(@(t,Y) NuRegTcellmodel(t,Y, params), tspan, init_cond, options);

% Relabel to easily keep track of desired compartments
AL = Y(:,1);
EL = Y(:,2);
RL = Y(:,3);

disp(T(end)./365);

% Plot EL vs RL  for low nu and the line RL = -EL + kR (steady state set)
figure(1)
subplot(1,2,1)
plot(EL./scale, RL./scale, 'LineWidth', 2, 'Color', 'b', 'LineStyle', '--')
plot([kE./scale 0], [0 kR./scale], '--', 'Color', [0.5, 0.5, 0.5], 'HandleVisibility', 'on', 'linewidth', 1)
% Plot specifications
xlabel('E_L \times 10^7', 'FontSize', 17)
ylabel('R_L \times 10^7', 'FontSize', 17)
legend('High \nu', 'Low \nu',...
'R_L = -E_L + k_R', 'FontSize', 17)
% xlim([0 0.025])
hold on
x = linspace(0, 3, 10); % Define x over a range, e.g., -10 to 10
y = x;                     % y = x
plot(x, y, 'k', 'LineWidth', 1.3, 'HandleVisibility', 'off'); % Plot the line in red with width 2

title('Late Disease Onset, R_L(0) = 10', 'FontSize', 17)

figure(2)
subplot(2,2,2)
rate1 = lambdaR.*AL;
rate2 = (omegaR).*(1 - (EL + RL)/C).*RL - phiR.*RL - deltaR.*RL;
plot(T./365, rate1./scale, 'LineWidth', 2, 'Color', 'b', 'LineStyle', '--')
hold on
plot(T./365, rate2./scale, 'LineWidth', 2, 'Color', 'r', 'LineStyle', '-')
% legend('First Rate for RL - with AL', 'Second Rate for RL - with carying capicity')
set(gca, 'Xscale', 'log')
title('low RL, low nu')
ylim([-5 10])



%% High RL
RL_initcond = 1e6;

% Save Initital Conditions in a Vector
init_cond = [AL_initcond EL_initcond RL_initcond AP_initcond ...
             EP_initcond RP_initcond B_initcond]';

% Specify Parameter Values
r1 = r1_vals(2);
r2 = r2_vals(2);
nu = nu_vals(1);

% Save parameters in a vector
params = [phiA deltaA lambdaEL omegaEL phiE deltaE lambdaR omegaR C ...
      phiR deltaR kappa r1 r2 nu]';

% Run the Model
options = odeset('Events', @(t, Y) events(t, Y, NuRegTcellmodel(t,Y, params), B_initcond));
[T,Y] = ode23s(@(t,Y) NuRegTcellmodel(t,Y, params), tspan, init_cond, options);

% Relabel to easily keep track of desired compartments
AL = Y(:,1);
EL = Y(:,2);
RL = Y(:,3);

disp(T(end)./365)

% Plot EL vs RL for high nu
figure(1)
subplot(1,2,2)
plot(EL./scale, RL./scale, 'LineWidth', 2, 'Color', 'b', 'LineStyle', '-')
xlabel('E_L \times 10^7', 'FontSize', 17)
ylabel('R_L \times 10^7', 'FontSize', 17)
% ylim([1.5 3])
hold on
x = linspace(0, 10, 3); % Define x over a range, e.g., -10 to 10
y = x;                     % y = x
plot(x, y, 'k', 'LineWidth', 1.3, 'HandleVisibility', 'off'); % Plot the line in red with width 2

figure(2)
subplot(2,2,3)
rate1 = lambdaR.*AL;
rate2 = (omegaR).*(1 - (EL + RL)/C).*RL - phiR.*RL - deltaR.*RL;
plot(T./365, rate1./scale, 'LineWidth', 2, 'Color', 'b', 'LineStyle', '--')
hold on
plot(T./365, rate2./scale, 'LineWidth', 2, 'Color', 'r', 'LineStyle', '-')
% legend('First Rate for RL - with AL', 'Second Rate for RL - with carying capicity')
set(gca, 'Xscale', 'log')
title('High RL, high nu')
ylim([-5 10])


% Plotting Low nu Value
nu = nu_vals(2);

% Save parameters in a vector
params = [phiA deltaA lambdaEL omegaEL phiE deltaE lambdaR omegaR C ...
      phiR deltaR kappa r1 r2 nu]';

% Run the Model
options = odeset('Events', @(t, Y) events(t, Y, NuRegTcellmodel(t,Y, params), B_initcond));
[T,Y] = ode23s(@(t,Y) NuRegTcellmodel(t,Y, params), tspan, init_cond, options);

% Relabel to easily keep track of desired compartments
AL = Y(:,1);
EL = Y(:,2);
RL = Y(:,3);

disp(T(end)./365);

% Plot EL vs RL  for low nu and the line RL = -EL + kR (steady state set)
figure(1)
subplot(1,2,2)
plot(EL./scale, RL./scale, 'LineWidth', 2, 'Color', 'b', 'LineStyle', '--')
plot([kE./scale 0], [0 kR./scale], '--', 'Color', [0.5, 0.5, 0.5], 'HandleVisibility', 'on', 'linewidth', 1)
% Plot specifications
xlabel('E_L \times 10^7', 'FontSize', 17)
ylabel('R_L \times 10^7', 'FontSize', 17)
legend('High \nu', 'Low \nu',...
'R_L = -E_L + k_R', 'FontSize', 17)
xlim([0 0.025])
title('Late Disease Onset, R_L(0) = 1 \times 10^6', 'FontSize', 17)

figure(2)
subplot(2,2,4)
rate1 = lambdaR.*AL;
rate2 = (omegaR).*(1 - (EL + RL)/C).*RL - phiR.*RL - deltaR.*RL;
plot(T./365, rate1./scale, 'LineWidth', 2, 'Color', 'b', 'LineStyle', '--')
hold on
plot(T./365, rate2./scale, 'LineWidth', 2, 'Color', 'r', 'LineStyle', '-')
legend('First Rate for RL - with AL', 'Second Rate for RL - with carying capicity')
set(gca, 'Xscale', 'log')
title('High RL, low nu')
ylim([-5 10])


set(gcf, 'Position', [100, 300, 1400, 500]);

% Save the figure as pdf
set(gcf, 'Units', 'Inches');
pos = get(gcf, 'Position');
set(gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)]);
% saveas(gcf, figurename); % Save Figure in Folder


%% Event Function
function [value, isterminal, direction] = events(~, y, ~, B_initcond)
    B = y(7);
    % disp(['B: ', num2str(B), ', Threshold: ', num2str(0.2 * B_initcond)]);
    value = B - 0.2*B_initcond;  % Condition: stop when dydt is close to zero
    isterminal = 1;  % Stop the integration
    direction = -1;   % Detect all zeros (both rising and falling)
end