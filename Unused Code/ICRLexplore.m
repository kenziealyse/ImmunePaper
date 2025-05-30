 
% CLEAR THE WORKSPACE
clc
clear
close all

% Set Model Conditions
tspan = 0:.1:70*365;

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

AL_initcond = 0;
EL_initcond = 10;
AP_initcond = 0;
EP_initcond = 0;
RP_initcond = 0;
B_initcond = 1*10^6;

N = 100;
RL_initconds = [10 10e6 10e10];%linspace(0, 1e6, N);


% Same or diff r1 values?
same = 1;
% Legend or no?
plotlegend = 1;

if same == 1
    % Same r1 and r2
    r1_vals = [0.2*10^(-5), 0.8*10^(-5)];
    r2_vals = [0.2*10^(-5), 1*10^(-5)];
    figurename = 'Figures/timetodisease_RLinitconds_samer1andr2.pdf';
else
    % Different r1 and r2
    r1_vals = [0.1*10^(-5), 0.8*10^(-5)];
    r2_vals = [0.1*10^(-5), 0.5*10^(-5)];
    figurename = 'Figures/timetodisease_RLinitconds_diffr1andr2.pdf';
end

nu =   1*10^(-5);%0.01;
colors = {'r', 'b', 'g'};


T_vals = zeros(length(RL_initconds), 1);


for j = 2:2

    r1 = r1_vals(j);
    r2 = r2_vals(j);

    params = [phiA deltaA lambdaEL omegaEL phiE deltaE lambdaR omegaR C ...
          phiR deltaR kappa r1 r2 nu]';

    for i = 1:length(RL_initconds)
    
        RL_initcond = RL_initconds(i)
              
    
        init_cond = [AL_initcond EL_initcond RL_initcond AP_initcond ...
             EP_initcond RP_initcond B_initcond]';
    
    
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
    
            T_vals(i) = T(end)./365

            figure(1)
    plot(T./365, (EL./1e7)./(RL./1e7), 'LineWidth', 1.3, 'Color', colors{i})
    hold on
    set(gca, 'Yscale', 'log')
    set(gca, 'Xscale', 'log')
    xlabel('Time, years', 'FontSize', 17)
    ylabel('EL/RL', 'FontSize', 17)
    legend('RL(0) = 10', 'RL(0) = 1e6', 'FontSize', 20)
    
    end

    

end

%% Low Nu
r1_vals = [0.2*10^(-5), 0.8*10^(-5)];
r2_vals = [0.2*10^(-5), 1*10^(-5)];

nu =   1*10^(-5);
colors = {'r', 'b'};

T_vals = zeros(length(RL_initconds), 1);


% for j = 2:2
% 
%     r1 = r1_vals(j);
%     r2 = r2_vals(j);
% 
%     params = [phiA deltaA lambdaEL omegaEL phiE deltaE lambdaR omegaR C ...
%           phiR deltaR kappa r1 r2 nu]';
% 
%     for i = 1:length(RL_initconds)
% 
%         RL_initcond = RL_initconds(i);
% 
% 
%         init_cond = [AL_initcond EL_initcond RL_initcond AP_initcond ...
%              EP_initcond RP_initcond B_initcond]';
% 
% 
%             % Run the Model
%             options = odeset('Events', @(t, Y) events(t, Y, NuRegTcellmodel(t,Y, params), B_initcond));
%             [T,Y] = ode23s(@(t,Y) NuRegTcellmodel(t,Y, params), tspan, init_cond, options);
% 
%             % Relabel to easily keep track of compartments
%             % LN
%             AL = Y(:,1);
%             EL = Y(:,2);
%             RL = Y(:,3);
% 
%             % Pancreas
%             AP = Y(:,4);
%             EP = Y(:,5);
%             RP = Y(:,6);
%             B = Y(:,7); 
% 
%             T_vals(i) = T(end)./365;
% 
%     end
% 
%     figure(1)
%     plot(T./365, RL./1e7, 'LineWidth', 1.3, 'Color', 'r', 'LineStyle', '--')
%     plot(T./365, EL./1e7, 'LineWidth', 1.3, 'Color', 'b', 'LineStyle', '--')
%     set(gca, 'Yscale', 'log')
%     set(gca, 'Xscale', 'log')
%     legend('high nu -RL', 'high nu EL', 'low nu -RL', 'low nu EL', 'FontSize', 20)
%     xlabel('Time')
% 
% end
% disp(T_vals)
% 
% if plotlegend == 1
%     legend('Early Disease Onset - high \nu', 'Late Disease Onset - high \nu',...
%         'Early Disease Onset - low \nu', 'Late Disease Onset - low \nu','FontSize', 17, 'Location', 'Best');
% else
%     %nothing
% end

% yline(70, 'Color', [0.5 0.5 0.5], 'LineStyle', '--', 'HandleVisibility', 'off', 'LineWidth', 1.2);
set(gcf, 'Position', [100, 300, 900, 600]);

% Save the figure as pdf
set(gcf, 'Units', 'Inches');
pos = get(gcf, 'Position');
set(gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)]);
saveas(gcf, figurename); % Save Figure in Folder


function [value, isterminal, direction] = events(~, y, ~, B_initcond)
    B = y(7);
    % disp(['B: ', num2str(B), ', Threshold: ', num2str(0.2 * B_initcond)]);
    value = B - 0.2*B_initcond;  % Condition: stop when dydt is close to zero
    isterminal = 1;  % Stop the integration
    direction = -1;   % Detect all zeros (both rising and falling)
end