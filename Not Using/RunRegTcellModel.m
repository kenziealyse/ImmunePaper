% CLEAR THE WORKSPACE
clc
clear
close all

% Set Model Conditions
tspan = 0:.1:1000;

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
sigma = params(12);
alpha = params(13);
deltaP = params(14);
kappa = params(15);
r1 = params(16);
r2 = params(17);

params = [phiA deltaA lambdaEL omegaEL phiE deltaE lambdaR omegaR C ...
          phiR deltaR sigma alpha deltaP kappa r1 r2]';

AL_initcond = 0;
EL_initcond = 1000;
RL_initcond = 0;
AP_initcond = 0;
P_initcond =  10;
EP_initcond = 0;
RP_initcond = 0;
B_initcond = 1*10^6;

init_cond = [AL_initcond EL_initcond RL_initcond AP_initcond P_initcond ...
             EP_initcond RP_initcond B_initcond]';

% Run the Model
options = odeset('Events', @(t, Y) events(t, Y, RegTcellModel(t,Y, params)));
[T,Y] = ode23s(@(t,Y) RegTcellModel(t,Y, params), tspan, init_cond, options);

% Relabel to easily keep track of compartments
% LN
AL = Y(:,1);
EL = Y(:,2);
RL = Y(:,3);

% Pancreas
AP = Y(:,4);
P = Y(:,5);
EP = Y(:,6);
RP = Y(:,7);
B = Y(:,8);

% Overlay
f = figure(1);
f.Position(3) = 1000;
f.Position(4) = 800;

T = T./365;

subplot(2,5,1)
plot(T, AL, 'LineWidth', 1.5, 'Color', 'b')
title('APC - LN')
xlabel('Time, years')
hold on

subplot(2,5,2)
plot(T, EL, 'LineWidth', 1.5, 'Color', 'b')
title('E - LN')
xlabel('Time, years')
hold on

subplot(2,5,3)
plot(T, RL, 'LineWidth', 1.5, 'Color', 'b')
title('Tregs - LN')
xlabel('Time, years')
hold on

subplot(2,5,4)
plot(T, AP, 'LineWidth', 1.5, 'Color', 'r')
title('APC - Pancreas')
xlabel('Time, years')
hold on

subplot(2,5,5)
plot(T, P, 'LineWidth', 1.5, 'Color', 'r')
title('Peptide - Pancreas')
xlabel('Time, years')
hold on

subplot(2,5,6)
plot(T, EP, 'LineWidth', 1.5, 'Color', 'r')
title('E - Pancreas')
xlabel('Time, years')
hold on

subplot(2,5,7)
plot(T, RP, 'LineWidth', 1.5, 'Color', 'r')
title('Tregs - Pancreas')
xlabel('Time, years')
hold on

subplot(2,5,8)
plot(T, B, 'LineWidth', 1.5, 'Color', 'r')
title('Beta Cells - Pancreas')
xlabel('Time, years')
hold on


function [value, isterminal, direction] = events(~, Y, dydt)
    value = dydt(8);  % Condition: stop when dydt is close to zero
    isterminal = 1;  % Stop the integration
    direction = 0;   % Detect all zeros (both rising and falling)
end
