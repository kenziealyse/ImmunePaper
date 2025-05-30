% CLEAR THE WORKSPACE
% clc
% clear
% close all

Init_B = 1*10^6;

params = LoadParameters();

phiA = params(1);
deltaA = [10*params(2):10:100*params(2)];
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

earlydiseasepoint = [0.2*10^(-5), 0.2*10^(-5)];
latediseasepoint = [0.8*10^(-5), 1*10^(-5)];

r1 = latediseasepoint(1);
r2 = latediseasepoint(2);
nu = 0.01;%0.25*10^(-5);%1*10^(-5);%

params = [phiA deltaA lambdaEL omegaEL phiE deltaE lambdaR omegaR C ...
                  phiR deltaR kappa r1 r2 nu]';

%% Pre Dose

AL_initcond = 0;
EL_initcond = 10;
RL_initcond = 1e6;
AP_initcond = 0;
EP_initcond = 0;
RP_initcond = 0;
B_initcond = 1*10^6;

init_cond = [AL_initcond EL_initcond RL_initcond AP_initcond ...
EP_initcond RP_initcond B_initcond]';

tspan1 = 0:0.01:150*365;
  

% Run the Model
options = odeset('Events', @(t, Y) events(t, Y, NuRegTcellmodel(t,Y, params), B_initcond));
[T,Y] = ode23s(@(t,Y) NuRegTcellmodel(t,Y, params), tspan1, init_cond, options);

% Relabel Compartments to Easily Keep Track
AL = Y(:,1);
EL = Y(:,2);
RL = Y(:,3);
AP = Y(:,4);
EP = Y(:,5);
RP = Y(:,6);
B = Y(:,7);
T(end)./365
% Plot APC term for high and low nu
APCrate = (lambdaEL.*AL)./(1 + r1.*RL);

% figure(2)
% subplot(1,2,1)
% plot(T./365, APCrate, 'Color', 'r')
% hold on
% set(gca, 'Xscale', 'log')
% title('APC activation Rate')
% 
% % Plot proliferation term for high and low nu
% prolifrate = omegaEL.*(1 - (EL + RL)./C).*EL;
% 
% subplot(1,2,2)
% plot(T./365, prolifrate, 'Color', 'r')
% hold on
% set(gca, 'Xscale', 'log')
% title('Proliferation Rate')

figure(2)
plot(T./365, AP, 'Color', 'b')
hold on
set(gca, 'Xscale', 'log')
title('APC in LN')

function [value, isterminal, direction] = events(~, y, ~, B_initcond)
    B = y(7);
%     disp(['B: ', num2str(B), ', Threshold: ', num2str(0.2 * B_initcond)]);
    value = B - .2*B_initcond;  % Condition: stop when dydt is close to zero
    isterminal = 1;  % Stop the integration
    direction = -1;   % Detect all zeros (both rising and falling)
end