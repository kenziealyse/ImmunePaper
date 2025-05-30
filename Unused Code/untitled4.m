 
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

% Plot EL vs RL  for low nu and the line RL = -EL + kR (steady state set)
figure(1)
plot([kE./scale 0], [0 kR./scale], '--', 'Color', [0.5, 0.5, 0.5], 'HandleVisibility', 'on', 'linewidth', 2)
hold on
x = linspace(0, 3, 100); % Define x over a range, e.g., -10 to 10
y = x;                     % y = x
plot(x, y, 'k', 'LineWidth', 1); % Plot the line in red with width 2
xlabel('E_L \times 10^7', 'FontSize', 17)
ylabel('R_L \times 10^7', 'FontSize', 17)
legend('R_L = -E_L + k_R', 'y=x', 'FontSize', 17)

%% Event Function
function [value, isterminal, direction] = events(~, y, ~, B_initcond)
    B = y(7);
    % disp(['B: ', num2str(B), ', Threshold: ', num2str(0.2 * B_initcond)]);
    value = B - 0.2*B_initcond;  % Condition: stop when dydt is close to zero
    isterminal = 1;  % Stop the integration
    direction = -1;   % Detect all zeros (both rising and falling)
end