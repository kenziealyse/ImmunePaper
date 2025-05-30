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
nu = 1e-5;

earlydiseasepoint = [0.2*10^(-5), 0.2*10^(-5)];
latediseasepoint = [0.8*10^(-5), 1*10^(-5)];

r1 = latediseasepoint(1);
r2 = latediseasepoint(2);

N = 100;
nu_vals = linspace(1*10^(-5), 2500*10^(-5), N); 

figurename = 'timetodiseasevsnu_ToyModel2_ND.pdf';

EL_initcond = 10/C;
RL_initcond = 10/C;
AP_initcond = 0;
B_initcond = (1e6.*nu.*lambdaR)/omegaR^2;

tspan1 = 0:0.01:(70*365)*omegaR;

init_cond = [EL_initcond RL_initcond AP_initcond B_initcond]';

params = [phiA deltaA lambdaEL omegaEL phiE deltaE lambdaR omegaR C ...
          phiR deltaR kappa r1 r2 nu]';

% Run the Model
options = odeset('Events', @(t, Y) Event(t, Y, B_initcond));
[T,Y] = ode23s(@(t,Y) ToyModel2_ND(t,Y, params), tspan1, init_cond, options);

% Save end time
T = T./(omegaR*365);

% Relabel Compartments to Easily Keep Track
EL = Y(:,1);
RL = Y(:,2);
AP = Y(:,3);
B = Y(:,4);

figure()
set(gca, 'Color', 'White')
subplot(2,2,1)
plot(T, EL.*C, 'LineWidth', 1.3, 'Color', 'k')
title('E', 'FontSize', 17)

subplot(2,2,2)
plot(T, RL.*C, 'LineWidth', 1.3, 'Color', 'k')
title('R', 'FontSize', 17)

subplot(2,2,3)
plot(T, (AP.*lambdaR)./(omegaR*C), 'LineWidth', 1.3, 'Color', 'k')
title('A', 'FontSize', 17)
set(gca, 'Yscale', 'log')

subplot(2,2,4)
plot(T, omegaR^2.*B./(nu.*lambdaR), 'LineWidth', 1.3, 'Color', 'k')
hold on
yline(0.2*omegaR^2.*B_initcond./(nu.*lambdaR), '--', 'LineWidth', 1.3, 'Color', [0.7 0.7 0.7])
title('\beta', 'FontSize', 17)
ylim([0 1e6])

function [value, isterminal, direction] = Event(~,Y, B_initcond)
    B = Y(4);
    value = B - 0.2*B_initcond;  % Condition: stop when dydt is close to zero
    isterminal = 1;  % Stop the integration
    direction = 0;   % Detect all zeros (both rising and falling)
end