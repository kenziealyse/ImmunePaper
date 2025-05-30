% CLEAR THE WORKSPACE
clc
clear
close all

% Set Model Conditions
tspan = 0:1:100000;

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

EL_initcondvals = 0:1000:10000;
RL_initcondvals = 0:10000:100000;

B_vals = zeros(length(RL_initcondvals), length(EL_initcondvals));

for i = 1:length(EL_initcondvals)

    for j = 1:length(RL_initcondvals)

        AL_initcond = 0;
        EL_initcond = EL_initcondvals(i);
        RL_initcond = RL_initcondvals(j);
        AP_initcond = 0;
        P_initcond =  0;
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

        B_vals(j,i) = B(end);

        subplot(1,2,1)
        plot(T./365, B, 'LineWidth', 1.5, 'Color', 'r')
        title('Beta Cells - Pancreas')
        xlabel('Time, years')
        hold on

    end  

end

% Define custom colormap
figure(1)
subplot(1,2,2)
imagesc(EL_initcondvals, RL_initcondvals, B_vals);
xlim([0 max(EL_initcondvals)])
ylim([0 max(RL_initcondvals)])
xlabel('Effector T cells in LN');
ylabel('Regulatory T cells in LN');
set(gca,'YDir','normal');  % Flip the y-axis to make it standardly oriented
title('Beta Cell Steady State Value');

% Define custom colormap
numColors = 100;
greenValues = linspace(0, 1, numColors);  % Green values range from 0 to 1
redValues = linspace(1, 0, numColors);    % Red values range from 1 to 0

% Create custom colormap without blue component
customColorMap = [redValues', greenValues', zeros(numColors, 1)];  % Red and green, blue is zero

% Apply custom colormap to the current figure
colormap(customColorMap);

% Show color bar with custom tick labels
colorbar;
caxis([0, max(B)]);  % Set color axis limits based on data range

function [value, isterminal, direction] = events(~, ~, dydt)
    B = dydt(8);
    value = B - 1e-2;  % Condition: stop when dydt is close to zero
    isterminal = 1;  % Stop the integration
    direction = 0;   % Detect all zeros (both rising and falling)
end