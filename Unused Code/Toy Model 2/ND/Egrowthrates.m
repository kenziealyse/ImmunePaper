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
nu = 0.01;

earlydiseasepoint = [0.2*10^(-5), 0.2*10^(-5)];
latediseasepoint = [0.8*10^(-5), 1*10^(-5)];

r1 = latediseasepoint(1);
r2 = latediseasepoint(2);

params = [phiA deltaA lambdaEL omegaEL phiE deltaE lambdaR omegaR C ...
              phiR deltaR kappa r1 r2 nu]';

N = 100;
nu_vals = [1e-5, 0.01];%linspace(1*10^(-5), 2500*10^(-5), N); 
colors = ['r', 'k'];
linestyles = {'-', '--'};
figurename = 'timetodiseasevsnu_ToyModel2_ND.pdf';

EL_initcond = 10/C;
RL_initcond = 1e6/C;
AP_initcond = 0;
B_vals = (1e6.*nu_vals.*lambdaR)/omegaR^2;

tspan1 = 0:0.01:(70*365)*omegaR;

figure(1)

for i = 1:length(B_vals)

    B_initcond = B_vals(i);
    init_cond = [EL_initcond RL_initcond AP_initcond B_initcond]';
       
    
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

    rate1 = (AP)./(1 + r1*C.*RL);
    subplot(1,2,1)
    plot(T, rate1, 'LineWidth', 1.3, 'Color', colors(i), 'LineStyle', linestyles(i))
    title(['rate1, RL(0) = ', num2str(RL_initcond.*C)], 'FontSize', 17)
    hold on

    subplot(1,2,2)
    rate2 = (1 - (EL + RL)).*EL;
    plot(T, rate2, 'LineWidth', 1.3, 'Color', colors(i), 'LineStyle', linestyles(i))
    title('rate2', 'FontSize', 17)
    hold on
end
set(gcf, 'Color', 'White')

legend(['\nu = ', num2str(nu_vals(1))], ['\nu = ', num2str(nu_vals(2))])

function [value, isterminal, direction] = Event(~,Y, B_initcond)
    B = Y(4);
    value = B - 0.2*B_initcond;  % Condition: stop when dydt is close to zero
    isterminal = 1;  % Stop the integration
    direction = 0;   % Detect all zeros (both rising and falling)
end