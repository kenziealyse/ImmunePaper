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

AL_initcond = 0;
EL_initcond = 10;
RL_initcond = 10;
AP_initcond = 0;
EP_initcond = 0;
RP_initcond = 0;
B_initcond = 1e6;

tspan1 = 0:0.01:70*365;


init_cond = [AL_initcond EL_initcond RL_initcond AP_initcond ...
EP_initcond RP_initcond B_initcond]';

param_names = {'\phi_A', '\delta_A', '\lambda_{E}', '\Omega_E', '\phi_E', '\delta_E',...
    '\lambda_R', '\Omega_R', 'C', '\phi_R', '\delta_R', '\kappa'};

orig_params = [phiA deltaA lambdaEL omegaEL phiE deltaE lambdaR omegaR C ...
                      phiR deltaR kappa r1 r2 nu]';

for j = 1:length(param_names)
        params = orig_params;
        SA_upperbound = params(j) + 0.1*params(j);
        SA_lowerbound = params(j) - 0.1*params(j);
        SA_param = [SA_lowerbound + (SA_upperbound - SA_lowerbound) * rand(1, 100), params(j)];
        SA_param = unique(SA_param);
        SA_param = sort(SA_param);
        time = zeros(length(SA_param), 1);
        
        for i = 1:length(SA_param)
            % Run the Model
            params(j) = SA_param(i);
            options = odeset('Events', @(t, Y) PercentBetaCellMassEvent(t, Y, NuRegTcellmodel(t,Y, params), B_initcond));
            [T,Y] = ode23s(@(t,Y) NuRegTcellmodel(t,Y, params), tspan1, init_cond, options);
            
            % Save end time
            time(i) = T(end)./365;
        end
        % Relabel Compartments to Easily Keep Track
        AL = Y(:,1);
        EL = Y(:,2);
        RL = Y(:,3);
        AP = Y(:,4);
        EP = Y(:,5);
        RP = Y(:,6);
        B = Y(:,7);
        
        subplot(3,4,j)
        plot((SA_param - orig_params(j))/orig_params(j), time, 'LineWidth', 1.3, 'Color', 'k')
        index = find(SA_param == orig_params(j));
        orig_time = time(index);
        hold on
        plot(0, orig_time, '*', 'LineWidth', 1.5, 'Color', 'r');
        xterm = strcat(string(param_names(j)), ' fold change');
        xlabel(xterm)
        ylabel('Time to 20% Beta Cell Mass, years')
        ylim([50 70])
end

screenSize = get(0, 'ScreenSize'); % [left bottom width height]
set(gcf, 'Position', screenSize);  % set figure to full screen
% Save the figure as pdf
figurename = 'paramSA.pdf';
set(gcf, 'Units', 'Inches');
pos = get(gcf, 'Position');
set(gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)]);
saveas(gcf, figurename); % Save Figure in Folder