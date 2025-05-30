% CLEAR THE WORKSPACE
clc
clear
close all

% Set Model Conditions
tspan = 0:1:70*365;

params = LoadParameters();

phiA = params(1);
deltaA = params(2);
lambdaEL = params(3);
omegaEL = params(4);
phiE = params(5);
deltaE = params(6);
lambdaR = params(7);
omegaR_vals = [params(8), 1.1*params(8), 0.9*params(8)];
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

indvplots = 0;
indvplots_compare = 1;
heatmapplots = 0;

figurename = 'lownuloweromega.pdf';

init_cond = [AL_initcond EL_initcond RL_initcond AP_initcond ...
         EP_initcond RP_initcond B_initcond]';


N = 300; % Desired number of points
r1_vals =  5.5*10^(-4);%linspace(0, 0.6*10^(-3), N); %[0.2*10^(-4), 2.2*10^(-4), 5.5*10^(-4)];%
nu_vals =   1*10^(-5);%linspace(2.5*10^(-6), 0.05, N); %[0.01];%
colors = {'r', 'g', 'r'};
linestyles = {'-', '-', '--'};

T_vals = zeros(length(nu_vals), length(r1_vals));
lambda = zeros(length(nu_vals), length(r1_vals));

for k = 1:length(omegaR_vals)
    omegaR = omegaR_vals(k)
for i = 1:length(r1_vals)

    r1 = r1_vals(i)

    

    for j = 1:length(nu_vals)
        
        nu = nu_vals(j);

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


        T_vals(j,i) = T(end)./365;
        disp(T(end)./365)

        lambda(j,i) = kappa*max(EP)/(1+r2*min(RP));

        if indvplots_compare == 1
            colors{k}
            figure(5)
            subplot(2,4,1)
            plot(T./365,AL, 'Color', colors{k}, 'LineWidth', 1.3, 'LineStyle', linestyles{k}) 
            hold on
            set(gca, 'Xscale', 'log')
            ylabel('A_L, cells', 'FontSize', 15)
            xlabel('Time, years', 'FontSize', 15)
            set(gca, 'Yscale', 'log')
            ylim([10e-1 10e8])
            xlim([0 10e2])
            yticks([10e-1 10e1 10e3 10e5 10e7])

            subplot(2,4,2)
            plot(T./365,EL, 'Color', colors{k}, 'LineWidth', 1.3, 'LineStyle', linestyles{k})  
            hold on
            set(gca, 'Xscale', 'log')
            ylabel('E_L, cells', 'FontSize', 15)
            xlabel('Time, years', 'FontSize', 15)
            legend('k_R = k_E', 'k_R > k_E', 'k_R < k_E', 'FontSize', 13)
            set(gca, 'Yscale', 'log')
            ylim([10e-1 10e8])
            xlim([0 10e2])
            yticks([10e-1 10e1 10e3 10e5 10e7])

            subplot(2,4,3)
            plot(T./365,RL, 'Color', colors{k}, 'LineWidth', 1.3, 'LineStyle', linestyles{k}) 
            hold on
            set(gca, 'Xscale', 'log')
            ylabel('R_L, cells', 'FontSize', 15)
            xlabel('Time, years', 'FontSize', 15)
            set(gca, 'Yscale', 'log')
            ylim([10e-1 10e8])
            xlim([0 10e2])
            yticks([10e-1 10e1 10e3 10e5 10e7])

            subplot(2,4,5)
            plot(T./365,AP, 'Color', colors{k}, 'LineWidth', 1.3, 'LineStyle', linestyles{k})
            hold on
            set(gca, 'Xscale', 'log')
            ylabel('A_P, cells', 'FontSize', 15)
            xlabel('Time, years', 'FontSize', 15)
            set(gca, 'Yscale', 'log')
            ylim([10e-1 10e8])
            xlim([0 10e2])
            yticks([10e-1 10e1 10e3 10e5 10e7])

            subplot(2,4,6)
            plot(T./365,EP, 'Color', colors{k}, 'LineWidth', 1.3, 'LineStyle', linestyles{k}) 
            hold on
            set(gca, 'Xscale', 'log')
            ylabel('E_P, cells', 'FontSize', 15)
            xlabel('Time, years', 'FontSize', 15)
            set(gca, 'Yscale', 'log')
            ylim([10e-1 10e8])
            xlim([0 10e2])
            yticks([10e-1 10e1 10e3 10e5 10e7])

            subplot(2,4,7)
            plot(T./365,RP, 'Color', colors{k}, 'LineWidth', 1.3, 'LineStyle', linestyles{k}) 
            hold on
            set(gca, 'Xscale', 'log')
            ylabel('R_P, cells', 'FontSize', 15)
            xlabel('Time, years', 'FontSize', 15)
            set(gca, 'Yscale', 'log')
            ylim([10e-1 10e8])
            xlim([0 10e2])
            yticks([10e-1 10e1 10e3 10e5 10e7])

            subplot(2,4,8)
            plot(T./365,B, 'Color', colors{k}, 'LineWidth', 1.3, 'LineStyle', linestyles{k}) 
            hold on
            yline(.2*B_initcond, '--', 'Color', 'k', 'LineWidth', 1.3);
            set(gca, 'Xscale', 'log')
            ylabel('B, cells', 'FontSize', 15)
            xlabel('Time, years', 'FontSize', 15)
            set(gca, 'Yscale', 'log')
            ylim([10^5 10e5])
            xlim([0 10])
            % Define custom y-ticks
            yticks([10^2, 10^4, 0.2 * B_initcond, 10^6]);
            % Define custom y-tick labels
            yticklabels({'10^2', '10^4', '20% B_0', '10^6'});
        end

    end  
 

end

end

set(gcf, 'Position', [100, 300, 1400, 500]);
    % set(gcf, 'Position', [100, 300, 800, 500]);

% Save the figure as pdf
set(gcf, 'Units', 'Inches');
pos = get(gcf, 'Position');
set(gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)]);
saveas(gcf, figurename); % Save Figure in Folder
   
function [value, isterminal, direction] = events(~, y, ~, B_initcond)
    B = y(7);
%     disp(['B: ', num2str(B), ', Threshold: ', num2str(0.2 * B_initcond)]);
    value = B - 0.2*B_initcond;  % Condition: stop when dydt is close to zero
    isterminal = 1;  % Stop the integration
    direction = -1;   % Detect all zeros (both rising and falling)
end