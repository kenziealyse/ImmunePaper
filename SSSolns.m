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
epsilon = params(9);
phiR = params(10);
deltaR = params(11);
sigma = params(12);
alpha = params(13);
deltaP = params(14);
kappa = params(15);
r1 = params(16);
r2 = params(17);

params = [phiA deltaA lambdaEL omegaEL phiE deltaE lambdaR omegaR epsilon ...
          phiR deltaR sigma alpha deltaP kappa r1 r2]';

EL_initcondvals = 100:10000:10000000;

for i = 1:length(EL_initcondvals)

    AL_initcond = 0;
    EL_initcond = 0;
    RL_initcond = EL_initcondvals(i);
    AP_initcond = 100;
    EP_initcond = 0;%EL_initcondvals(i);
    RP_initcond = 0;
    B_initcond = 1*10^6;
    
    init_cond = [AL_initcond EL_initcond RL_initcond AP_initcond ...
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
    EP = Y(:,5);
    RP = Y(:,6);
    B = Y(:,7);
    
    % Overlay
    f = figure(1);
    f.Position(3) = 1500;
    f.Position(4) = 600;
    
    T = T./365;
    
    if B(end) >= .2*B_initcond
        B_vals(i) = 1;
    else
        B_vals(i) = 0;
    end
    
    subplot(1,7,1)
    plot(T, AL, 'LineWidth', 1.5, 'Color', 'b')
    title('APC - LN')
    xlabel('Time, years')
    hold on
    
    subplot(1,7,2)
    plot(T, EL, 'LineWidth', 1.5, 'Color', 'b')
    title('E - LN')
    xlabel('Time, years')
    hold on
    
    subplot(1,7,3)
    plot(T, RL, 'LineWidth', 1.5, 'Color', 'b')
    title('Tregs - LN')
    xlabel('Time, years')
    hold on
    
    subplot(1,7,4)
    plot(T, AP, 'LineWidth', 1.5, 'Color', 'r')
    title('APC - Pancreas')
    xlabel('Time, years')
    hold on
    
    subplot(1,7,5)
    plot(T, EP, 'LineWidth', 1.5, 'Color', 'r')
    title('E - Pancreas')
    xlabel('Time, years')
    hold on
    
    subplot(1,7,6)
    plot(T, RP, 'LineWidth', 1.5, 'Color', 'r')
    title('Tregs - Pancreas')
    xlabel('Time, years')
    hold on
    
    subplot(1,7,7)
    plot(T, B./B_initcond, 'LineWidth', 1.5, 'Color', 'r')
    title('Beta Cells - Pancreas')
    xlabel('Time, years')
    ylim([0 1])
    hold on

    figure(3)
    plot(EL(end)./10e7, RL(end)./10e7, 'o', 'LineWidth', 1.5, 'Color', 'r')
    hold on
%     ylim([0 1])
    xlabel('EL')
    ylabel('RL')

end

figure(2)
histogram(B_vals,  [-0.25, 0.25, 0.75, 1.35], 'FaceColor', 'k', 'EdgeColor', 'black')
xticks([0 1]);  % Set the tick positions
xticklabels({'Disease', 'Healthy'});  % Set the tick labels

function [value, isterminal, direction] = events(~, ~, dydt)
    B = dydt(7);
    value = B - 1e-2;  % Condition: stop when dydt is close to zero
    isterminal = 1;  % Stop the integration
    direction = 0;   % Detect all zeros (both rising and falling)
end