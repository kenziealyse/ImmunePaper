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

N = 100;
nu_vals = linspace(1*10^(-5), 2500*10^(-5), N); 
figurename = 'timetodiseasevsnu_smallmodel.pdf';

EL_initcond = 10/C;
RL_initcond = 10/C;
AP_initcond = 0;
B_vals = (1e6.*nu_vals.*lambdaR)/omegaR^2;

tspan1 = 0:0.01:(70*365)*omegaR;

% Preallocate Space
time = zeros(length(B_vals), 1);
  
for i = 1:length(B_vals)

    B_initcond = B_vals(i);  

    init_cond = [EL_initcond RL_initcond AP_initcond B_initcond]';


    params = [phiA deltaA lambdaEL omegaEL phiE deltaE lambdaR omegaR C ...
                  phiR deltaR kappa r1 r2]';

    % Run the Model
    options = odeset('Events', @(t, Y) RatioEvent(t, Y, B_initcond));
    [T,Y] = ode23s(@(t,Y) SmallModel2_ND(t,Y, params), tspan1, init_cond, options);

    % Save end time
    time(i) = T(end)./(365);
    

end
figure(1)
plot(nu_vals, time./omegaR, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-')
hold on

%% Increase RL(0)
RL_initcond = 1e6/C;

% Preallocate Space
time = zeros(length(B_vals), 1);
  
for i = 1:length(B_vals)

    B_initcond = B_vals(i);
    init_cond = [EL_initcond RL_initcond AP_initcond B_initcond]';

    params = [phiA deltaA lambdaEL omegaEL phiE deltaE lambdaR omegaR C ...
                  phiR deltaR kappa r1 r2]';

     % Run the Model
    options = odeset('Events', @(t, Y) RatioEvent(t, Y, B_initcond));
    [T,Y] = ode23s(@(t,Y) SmallModel2_ND(t,Y, params), tspan1, init_cond, options);

    % Save end time
    time(i) = T(end)./(365);

end

figure(1)
plot(nu_vals, time./omegaR, 'Color', 'r', 'LineWidth', 1.5, 'LineStyle', '--')
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

function [value, isterminal, direction] = RatioEvent(~,Y, B_initcond)
    % ratio = Y(5);
    % value = ratio - log(5)/kappa;  % Condition: stop when dydt is close to zero
    B = Y(4);
    value = B - 0.2.*B_initcond;
    isterminal = 1;  % Stop the integration
    direction = 0;   % Detect all zeros (both rising and falling)
end