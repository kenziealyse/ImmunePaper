% Clear the workspace
clc 
clear
close all
% Initialize the counter as a global variable
global event_count;
event_count = 0;
% Load Parameter Values
params = LoadParameters();

% Specify Parameter Values
phiA = params(1);
deltaA_vals = [0:.1:100*params(2)];
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

% Specify Initial Conditions
AL_initcond = 0;
EL_initcond = 10;
RL_initcond = 1e6;
AP_initcond = 0;
EP_initcond = 0;
RP_initcond = 0;
B_initcond = 1*10^6;

tspan = 0:.1:70*365;

kR = 1/((omegaR)/(C*(omegaR - (phiR+deltaR))));
kE = 1/((omegaEL)/(C*(omegaEL - (phiE+deltaE))));
scale = 1e7;

nu_vals = 1*10^(-5);%sort([1*10^(-5):0.0001:0.025]);


% Save Initital Conditions in a Vector
init_cond = [AL_initcond EL_initcond RL_initcond AP_initcond ...
             EP_initcond RP_initcond B_initcond]';
j = 0;

for i = 1:length(deltaA_vals)
    
    nu = nu_vals;
    deltaA = deltaA_vals(i);
    % Save parameters in a vector
    params = [phiA deltaA lambdaEL omegaEL phiE deltaE lambdaR omegaR C ...
      phiR deltaR kappa r1 r2 nu]';
    % Run the Model
    options = odeset('Events', @(t, Y) events(t, Y, NuRegTcellmodel(t,Y, params), B_initcond));
    [T, Y, te, ye, ie]= ode23s(@(t,Y) NuRegTcellmodel(t,Y, params), tspan, init_cond, options);

    if length(te) > 1
        j = j+1;
    end
    % Relabel Compartments to Easily Keep Track
    AL = Y(:,1);
    EL = Y(:,2);
    RL = Y(:,3);
    AP = Y(:,4);
    EP = Y(:,5);
    RP = Y(:,6);
    B = Y(:,7);

    % figure(1)
    % subplot(1,2,1)
    % plot(T./365, EL./1e7, 'LineWidth', 1.3);
    % hold on
    % set(gca, 'Xscale', 'log')
    % xlabel('Time, years', 'FontSize', 17)
    % ylabel('E_L', 'FontSize', 17)
    % subplot(1,2,2)
    % plot(T./365, RL./1e7, 'LineWidth', 1.3);
    % xlabel('Time, years', 'FontSize', 17)
    % ylabel('R_L', 'FontSize', 17)
    % set(gca, 'Xscale', 'log')
    % % set(gca, 'Yscale', 'log')
    % legendEntries{i} = ['\nu = ', num2str(nu_vals(i))];
    % hold on
    
    EL_SS(i) = EL(end);
    T_SS(i) = T(end)./365;
    maxes(i) = max(EL);

    % figure(5)
    % plot(EL./1e7, RL./1e7, 'LineWidth', 1.3)
    % hold on
    % plot([kE./scale 0], [0 kR./scale], '--', 'Color', [0.5, 0.5, 0.5], 'HandleVisibility', 'on', 'linewidth', 1)
    % kE = 1/((omegaEL)/(C*(omegaEL - (phiE+deltaE))));
    % 
   % figure(6)
   %  rate = (lambdaEL.*AL)./(1+r1.*RL);
   %  plot(T./365, rate, 'LineWidth', 1.3);
   %  hold on
   %  rate2 = (lambdaR.*AL);
   %  plot(T./365, rate2, 'LineWidth', 1.3);
   %  hold on
   %  set(gca, 'Xscale', 'log')
   %  xlabel('Time, years', 'FontSize', 17)
   %  ylabel('Rates', 'FontSize', 17)
end

% figure(1)
% legend(legendEntries, 'FontSize', 20);

% figure(2)
% subplot(2,2,1)
% plot(nu_vals, EL_SS, 'LineWidth', 1.3)
% xlabel('\nu', 'FontSize', 17)
% ylabel('Steady State of EL', 'FontSize', 17)
% subplot(2,2,2)
% plot(nu_vals, T_SS, 'LineWidth', 1.3)
% xlabel('\nu', 'FontSize', 17)
% ylabel('Time to Disease, years', 'FontSize', 17)
% subplot(2,2,3)
% plot(nu_vals, maxes, 'LineWidth', 1.3)
% xlabel('\nu', 'FontSize', 17)
% ylabel('Max of EL', 'FontSize', 17)

figure(3)
plot(deltaA_vals, T_SS, 'LineWidth', 1.3, 'Color', 'k', 'LineStyle', '--')
hold on
% scatter(1*10^(-5), 34.6806, 'filled', 'MarkerEdgeColor', 'r', ...
%         'MarkerFaceColor', 'r', 'LineWidth',1.5)
% scatter(0.01, 64.0012, 'filled', 'MarkerEdgeColor', 'r', ...
%         'MarkerFaceColor', 'r', 'LineWidth',1.5)
xlabel('APC death rate', 'FontSize', 17)
ylabel('Time to Disease, years', 'FontSize', 17)


%% Function

%% Event Function
function [value, isterminal, direction] = events(~, Y, dydt, IC)

    EL = Y(2);
    RL = Y(3);
    B = Y(7);
% B = y(3);
    value = B - 0.2*IC;
    % disp(['B: ', num2str(B), ', Threshold: ', num2str(0.2 * B_initcond)]);
    % value = 0;%dydt(2);  % Condition: stop when dydt is close to zero
    % value = -(EL + RL) + kR;
    isterminal = 1;  % Stop the integration
    direction = -1;   % Detect all zeros (both rising and falling)


    
end