% CLEAR THE WORKSPACE
clc
clear
close all

params = LoadParameters();
phiA = params(1);
delta_A_orig = params(2);
deltaA_vals = 0:1:100*delta_A_orig;
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
nu = 1e-5;

figurename =['Figures\timetodiseasevsAPCdeathrate', num2str(nu),'.pdf'];

AL_initcond = 0;
EL_initcond = 10;
RL_initcond = 10;
AP_initcond = 0;
EP_initcond = 0;
RP_initcond = 0;
B_initcond = 1*10^6;

init_cond = [AL_initcond EL_initcond RL_initcond AP_initcond ...
EP_initcond RP_initcond B_initcond]';

tspan1 = 0:0.01:70*365;

% Preallocate Space
time = zeros(length(deltaA_vals), 1);
  
for i = 1:length(deltaA_vals)

    deltaA = deltaA_vals(i);

    params = [phiA deltaA lambdaEL omegaEL phiE deltaE lambdaR omegaR C ...
                  phiR deltaR kappa r1 r2 nu]';

    % Run the Model
    options = odeset('Events', @(t, Y) PercentBetaCellMassEvent(t, Y, NuRegTcellmodel(t,Y, params), B_initcond));
    [T,Y] = ode23s(@(t,Y) NuRegTcellmodel(t,Y, params), tspan1, init_cond, options);

    % Save end time
    time(i) = T(end)./365;

end

figure(1)
plot(deltaA_vals, time, 'Color', 'k', 'LineWidth', 1.3, 'LineStyle', '--')
hold on
ylabel('Time to disease', 'FontSize', 17)
xlabel('APC death rate, \delta_A', 'FontSize', 17)


%% Increase RL(0)
RL_initcond = 1e6;

init_cond = [AL_initcond EL_initcond RL_initcond AP_initcond ...
EP_initcond RP_initcond B_initcond]';

% Preallocate Space
time = zeros(length(deltaA_vals), 1);
  
for i = 1:length(deltaA_vals)

    deltaA = deltaA_vals(i);

    params = [phiA deltaA lambdaEL omegaEL phiE deltaE lambdaR omegaR C ...
                  phiR deltaR kappa r1 r2 nu]';

    % Run the Model
    options = odeset('Events', @(t, Y) PercentBetaCellMassEvent(t, Y, NuRegTcellmodel(t,Y, params), B_initcond));
    [T,Y] = ode23s(@(t,Y) NuRegTcellmodel(t,Y, params), tspan1, init_cond, options);

    % Save end time
    time(i) = T(end)./365;

end

figure(1)
plot(deltaA_vals, time, 'Color', 'r', 'LineWidth', 1.3, 'LineStyle', '--')
ylabel('Time to 20% Beta Cell Mass', 'FontSize', 17)
xlabel('APC death rate, \delta_A', 'FontSize', 17)
legend('R_L(0) = 10', 'R_L(0) = 1 \times 10^6', 'FontSize', 15, 'Location', 'best')


set(gcf, 'Position', [100, 300, 800, 500]);
set(gcf, 'Color', 'White')
% Save the figure as pdf
set(gcf, 'Units', 'Inches');
pos = get(gcf, 'Position');
set(gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)]);
saveas(gcf, figurename); % Save Figure in Folder