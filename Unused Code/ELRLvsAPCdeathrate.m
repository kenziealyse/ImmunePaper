% CLEAR THE WORKSPACE
clc
clear
% close all

params = LoadParameters();
phiA = params(1);
delta_A_orig = params(2);
deltaA_vals = delta_A_orig;%[3, 20];%0:1:100*delta_A_orig;
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
scale = 1e10;


earlydiseasepoint = [0.2*10^(-5), 0.2*10^(-5)];
latediseasepoint = [0.8*10^(-5), 1*10^(-5)];

r1 = latediseasepoint(1);
r2 = latediseasepoint(2);
nu = 0.01;%1*10^(-5);%

figurename =['Figures\timetodiseasevsAPCdeathrate', num2str(nu),'.pdf'];

AL_initcond = 0;
EL_initcond = 10;
RL_initcond = 10; 
AP_initcond = 0;
EP_initcond = 0;
RP_initcond = 0;
B_initcond = 1*10^6;

init_cond = [AL_initcond EL_initcond RL_initcond AP_initcond ...
EP_initcond RP_initcond B_initcond 0]';

tspan1 = 0:0.01:70*365;

% Preallocate Space
RLEL_ratio = zeros(length(deltaA_vals), 1);

kR = 1/((omegaR)/(C*(omegaR - (phiR+deltaR))));
kE = 1/((omegaEL)/(C*(omegaEL - (phiE+deltaE))));

  
for i = 1:length(deltaA_vals)

    deltaA = deltaA_vals(i);

    params = [phiA deltaA lambdaEL omegaEL phiE deltaE lambdaR omegaR C ...
                  phiR deltaR kappa r1 r2 nu]';

    % Run the Model
    options = odeset('Events', @(t, Y) PercentBetaCellMassEvent(t, Y, NuRegTcellmodel(t,Y, params), B_initcond));
    [T,Y] = ode23s(@(t,Y) NuRegTcellmodel(t,Y, params), tspan1, init_cond, options);
    
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
    ratio = Y(:,8);


    figure(2)
    plot(T./365, ratio, 'LineWidth', 1.3, 'Color', 'k', 'LineStyle', '--')
    hold on
    % xlim([0 0.00005])
    T(end)./365
    % set(gca, 'Xscale', 'log')

end

% legend('\delta_A = 3', '\delta_A = 20')
% plot([kE./scale 0], [0 kR./scale], '--', 'Color', [0.5, 0.5, 0.5], 'HandleVisibility', 'off', 'linewidth', 2)

% min_ratio = min(RLEL_ratio);
% max_ratio = max(RLEL_ratio);
% 
% figure(1)
% plot(deltaA_vals, RLEL_ratio./365, 'Color', 'k', 'LineWidth', 1.3, 'LineStyle', '--')
% hold on


%% Increase RL(0)
RL_initcond = 1e6;

init_cond = [AL_initcond EL_initcond RL_initcond AP_initcond ...
EP_initcond RP_initcond B_initcond 0]';

% Preallocate Space
RLEL_ratio = zeros(length(deltaA_vals), 1);

for i = 1:length(deltaA_vals)

    deltaA = deltaA_vals(i);

    params = [phiA deltaA lambdaEL omegaEL phiE deltaE lambdaR omegaR C ...
                  phiR deltaR kappa r1 r2 nu]';

    % Run the Model
    options = odeset('Events', @(t, Y) PercentBetaCellMassEvent(t, Y, NuRegTcellmodel(t,Y, params), B_initcond));
    [T,Y] = ode23s(@(t,Y) NuRegTcellmodel(t,Y, params), tspan1, init_cond, options);

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
    ratio = Y(:,8);

    % Save end RL:EL ratio
    [max_value, max_index] = max(AL);
    RLEL_ratio(i) = T(max_index);

    figure(2)
    plot(T./365, ratio, 'LineWidth', 1.3, 'Color', 'r', 'LineStyle', '--')
    hold on
    % xlim([0 0.00005])
    T(end)./365
    set(gca, 'Xscale', 'log')
    yline(log(5)/kappa, 'Color', 'k')

end

% legend('Low RL(0) - low nu', 'High RL(0) - low nu', 'FontSize', 17)
% 
% figure(1)
% plot(deltaA_vals, RLEL_ratio./365, 'Color', 'r', 'LineWidth', 1.3, 'LineStyle', '--')
% ylabel('Time of max(AL)', 'FontSize', 20)
% xlabel('APC death rate, \delta_A', 'FontSize', 20)
% legend('R_L(0) = 10', 'R_L(0) = 1 \times 10^6', 'FontSize', 20, 'Location', 'best')
% ax = gca;
% ax.FontSize = 20;
% 
% set(gcf, 'Position', [100, 300, 800, 500]);
% 
% % Save the figure as pdf
% set(gcf, 'Units', 'Inches');
% pos = get(gcf, 'Position');
% set(gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)]);
% saveas(gcf, figurename); % Save Figure in Folder