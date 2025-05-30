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

earlydiseasepoint = [0.2*10^(-5), 0.2*10^(-5)];
latediseasepoint = [0.8*10^(-5), 1*10^(-5)];

r1 = latediseasepoint(1);
r2 = latediseasepoint(2);

N = 10;
nu_vals =  linspace(1*10^(-5), 100, N); 

figurename = 'timetodiseasevsnu_smallmodel.pdf';

AL_initcond = 0;
EL_initcond = 10;
RL_initcond = 10;
AP_initcond = 0;
EP_initcond = 0;
RP_initcond = 0;
B_initcond = 1*10^6;

init_cond = [EL_initcond RL_initcond AP_initcond B_initcond]';

tspan1 = 0:0.1:70*365;
scale = 1e7;
% Preallocate Space
time = zeros(length(nu_vals), 1);

% Figure
f = figure(1);
f.Position(4) = 600;
f.Position(3) = 1100;
  
for i = 1:length(nu_vals)

    nu = nu_vals(i);

    params = [phiA deltaA lambdaEL omegaEL phiE deltaE lambdaR omegaR C ...
                  phiR deltaR kappa r1 r2 nu]';

    % Run the Model
    options = odeset('Events', @(t, Y) RatioEvent(t, Y, params, nu_vals(i)));
    [T,Y] = ode23s(@(t,Y) SmallModel2_ND(t,Y, params), tspan1, init_cond, options);
    % Relabel Compartments to Easily Keep Track
    EL = Y(:,1);
    RL = Y(:,2);
    AP = Y(:,3);
    B = Y(:,4);

    % figure(1)
    % subplot(1,2,1)
    % plot(T./365, EL, 'Color', 'b')
    % set(gca, 'Yscale', 'log')
    % title('E_L')
    % hold on
    % subplot(1,2,2)
    % plot(T./365, RL, 'Color', 'b')
    % set(gca, 'Yscale', 'log')
    % title('R_L')
    % hold on



    % Save end time
    time(i) = T(end)./365;
        [max_value, max_index] = max(AP);

    maxes(i) = T(max_index);


end
figure(5)
plot(nu_vals, maxes, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-')
hold on

figure(2)
plot(nu_vals, time, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-')
hold on

%% Increase RL(0)
 RL_initcond = 1e6;

init_cond = [EL_initcond RL_initcond AP_initcond B_initcond]';

% Preallocate Space
time = zeros(length(nu_vals), 1);
  
for i = 1:length(nu_vals)

    nu = nu_vals(i);

    params = [phiA deltaA lambdaEL omegaEL phiE deltaE lambdaR omegaR C ...
                  phiR deltaR kappa r1 r2 nu]';

     % Run the Model
    options = odeset('Events', @(t, Y) RatioEvent(t, Y, params, nu_vals(i)));
    [T,Y] = ode23s(@(t,Y) SmallModel2_ND(t,Y, params), tspan1, init_cond, options);

    % Relabel Compartments to Easily Keep Track
    EL = Y(:,1);
    RL = Y(:,2);
    AP = Y(:,3);
    B = Y(:,4);

    % Save end time
    time(i) = T(end)./365;

    rate1 = AP./(1 + r1.*C.*RL);
    rate2 = (1 - (EL + RL)).*EL;
figure(10)
subplot(1,2,1)
plot(T./365, nu.*EL.*B , 'Color', 'm', 'LineWidth', 1.2)
set(gca, 'Yscale', 'log')
hold on
subplot(1,2,2)
plot(T./365, AP,'Color', 'm', 'LineWidth', 1.2)
hold on
    figure(1)
    subplot(1,4,1)
    plot(T./365, EL, 'Color', 'r', 'LineWidth', 1.2)
    hold on
    % plot(T./365, rate2, 'Color', 'b', 'LineWidth', 1.2)
    title('E_L')
    hold on
    subplot(1,4,2)
    plot(T./365, RL, 'Color', 'r', 'LineWidth', 1.2)
    title('R_L')
    hold on
    subplot(1,4,3)
    plot(T./365, B./B(1))
    hold on
    yline(0.2*B(1)./B(1), '--', 'Color', 'k', 'LineWidth', 1.2)
    ylim([0 1])
    subplot(1,4,4)
    plot(T./365, AP, 'Color', 'r', 'LineWidth', 1.2)
    hold on
    [max_value, max_index] = max(AP);

    maxes(i) = T(max_index);
end
figure(5)
plot(nu_vals, maxes, 'Color', 'r', 'LineWidth', 1.5, 'LineStyle', '--')

figure(2)
plot(nu_vals, time, 'Color', 'r', 'LineWidth', 1.5, 'LineStyle', '--')
hold on
% xline(0.01)
ax = gca;
ax.FontSize = 21;
ylabel('Time to 20% Beta Cell Mass')
xlabel('\nu', 'FontSize', 30)
legend('R_L(0) = 10', 'R_L(0) = 1 \times 10^6', 'Location', 'best')

set(gcf, 'Position', [100, 300, 800, 500]);

% Save the figure as pdf
set(gcf, 'Units', 'Inches');
pos = get(gcf, 'Position');
set(gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)]);
saveas(gcf, figurename); % Save Figure in Folder

function [value, isterminal, direction] = RatioEvent(~,Y, params, nu)
    B = Y(4);
    value = B - 0.2.*1*10^6;

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

earlydiseasepoint = [0.2*10^(-5), 0.2*10^(-5)];
latediseasepoint = [0.8*10^(-5), 1*10^(-5)];

r1 = latediseasepoint(1);
r2 = latediseasepoint(2);

    % Relabel Compartments to Easily Keep Track
EL = Y(1);
RL = Y(2);
AP = Y(3);
B = Y(4);
    % value = AP/(1 + r1*C*RL) + (1 - (EL + RL))*EL;
    isterminal = 1;  % Stop the integration
    direction = 0;   % Detect all zeros (both rising and falling)
end