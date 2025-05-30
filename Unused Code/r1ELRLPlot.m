% CLEAR THE WORKSPACE
clc
clear
close all

% Set Model Conditions
tspan = 0:.01:150*365;

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
r2 = params(17);

AL_initcond = 0;
EL_initcond = 10;
RL_initcond = 10;
AP_initcond = 0;
EP_initcond = 0;
RP_initcond = 0;
B_initcond = 1*10^6;

r1_vals =   [0.2*10^(-4), 2.2*10^(-4), 5.5*10^(-4)];
nu =  0.01;

figurename = ['Figures/ELvsRLplotr1nuis', num2str(nu), '.pdf'];

scale = 1e7;

kR = 1/((omegaR)/(C*(omegaR - (phiR+deltaR))));
kE = 1/((omegaEL)/(C*(omegaEL - (phiE+deltaE))));

T_vals = zeros(length(r1_vals), 1);
lambda = zeros(length(r1_vals), 1);
init_cond = [AL_initcond EL_initcond RL_initcond AP_initcond ...
         EP_initcond RP_initcond B_initcond]';

colors = {'r', 'b', 'g', 'k'};


for i = 1:length(r1_vals)

    if i~= 4
        r1 = r1_vals(i);
    end
        params = [phiA deltaA lambdaEL omegaEL phiE deltaE lambdaR omegaR C ...
          phiR deltaR kappa r1 r2 nu]';
    
        % Run the Model
        options = odeset('Events', @(t, Y) events(t, Y, NuRegTcellmodel(t,Y, params), B_initcond));
        [T,Y] = ode23s(@(t,Y) NuRegTcellmodel(t,Y, params), tspan, init_cond, options);
        
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

        T_vals(i) = T(end)./365;

        disp(T(end)./365)

        figure(1)
        plot(EL./scale, RL./scale, 'LineWidth', 1.3, 'Color', colors{i})
        ylabel('R_L \times 10^7', 'FontSize', 17)
        hold on
        % scatter(EL(end)./scale,RL(end)./scale, 'Color', colors{i}, 'LineWidth', 5, 'Marker', '*', 'HandleVisibility', 'off') 
        

end

% disp(T_vals)

figure(1)
plot([kE./scale 0], [0 kR./scale], '--', 'Color', [0.5, 0.5, 0.5], 'HandleVisibility', 'on', 'linewidth', 2)
xlim([0 .1])
ylim([0 12])
legend('Early Disease Onset', 'Late Disease Onset', 'No Disease', 'R_L = -E_L + k_R', 'FontSize', 13)
        xlabel('E_L \times 10^7', 'FontSize', 17)
set(gcf, 'Position', [100, 300, 900, 600]);

% Save the figure as pdf
set(gcf, 'Units', 'Inches');
pos = get(gcf, 'Position');
set(gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)]);
saveas(gcf, figurename); % Save Figure in Folder


%% Low Nu
nu =  1*10^(-5);
r1 = 5.5*10^(-4);

figurename = ['Figures/ELvsRLplotr1nuis', num2str(nu), '.pdf'];

 params = [phiA deltaA lambdaEL omegaEL phiE deltaE lambdaR omegaR C ...
          phiR deltaR kappa r1 r2 nu]';
    
% Run the Model
options = odeset('Events', @(t, Y) events(t, Y, NuRegTcellmodel(t,Y, params), B_initcond));
[T,Y] = ode23s(@(t,Y) NuRegTcellmodel(t,Y, params), tspan, init_cond, options);

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

disp(T(end)./365)

figure(2)
plot(EL./scale, RL./scale, 'LineWidth', 1.3, 'Color', 'k')
ylabel('R_L \times 10^7', 'FontSize', 17)
hold on
% scatter(EL(end)./scale,RL(end)./scale, 'Color', 'k', 'LineWidth', 5, 'Marker', 'o', 'HandleVisibility', 'off') 
plot([kE./scale 0], [0 kR./scale], '--', 'Color', [0.5, 0.5, 0.5], 'HandleVisibility', 'on', 'linewidth', 2)
xlim([0 1])
ylim([0 4])
legend('Early Disease Onset', 'R_L = -E_L + k_R', 'FontSize', 13)
        xlabel('E_L \times 10^7', 'FontSize', 17)
set(gcf, 'Position', [100, 300, 900, 600]);

% Save the figure as pdf
set(gcf, 'Units', 'Inches');
pos = get(gcf, 'Position');
set(gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)]);
saveas(gcf, figurename); % Save Figure in Folder



function [value, isterminal, direction] = events(~, y, ~, B_initcond)
    B = y(7);
    % disp(['B: ', num2str(B), ', Threshold: ', num2str(0.2 * B_initcond)]);
    value = B - 0*B_initcond;  % Condition: stop when dydt is close to zero
    isterminal = 1;  % Stop the integration
    direction = -1;   % Detect all zeros (both rising and falling)
end