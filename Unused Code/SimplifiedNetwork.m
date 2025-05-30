% Clear the workspace
clc 
clear
close all

EL_initcond = 10;
AP_initcond = 0;
B_initcond = 1e6;

tspan = 0:.1:70*365;

nu = .25*1*10^(-5):.01:0.025;

init_cond = [EL_initcond, AP_initcond, B_initcond];


for i = 1:length(nu)
    % Run the Model
    options = odeset('Events', @(t, Y) events(t, Y, SimplifiedEqn(t,Y, nu(i)), B_initcond));
    [T,Y] = ode23s(@(t,Y) SimplifiedEqn(t,Y, nu(i)), tspan, init_cond, options);

    EL = Y(:,1);
    AP = Y(:,2);
    B = Y(:,3);

    figure(1)
    subplot(1,3,i)
    plot(T./365, EL./1e7, 'LineWidth', 1.3);
    set(gca, 'Xscale', 'log')
    % set(gca, 'Yscale', 'log')
    legendEntries{i} = ['\nu = ', num2str(nu(i))];
    hold on
    
    EL_SS(i) = EL(end);
    T_SS(i) = T(end)./365;
    maxes(i) = max(EL);
    
end

legend(legendEntries, 'FontSize', 20);

figure(2)
subplot(2,2,1)
plot(nu, EL_SS, 'LineWidth', 1.3)
xlabel('\nu', 'FontSize', 17)
ylabel('Steady State of EL', 'FontSize', 17)
subplot(2,2,2)
plot(nu, T_SS, 'LineWidth', 1.3)
xlabel('\nu', 'FontSize', 17)
ylabel('Time to Steady State, years', 'FontSize', 17)
subplot(2,2,3)
plot(nu, maxes, 'LineWidth', 1.3)
xlabel('\nu', 'FontSize', 17)
ylabel('Max of EL', 'FontSize', 17)

figure(3)
plot(nu, T_SS, 'LineWidth', 1.3)
xlabel('\nu', 'FontSize', 17)
ylabel('Time to Disease, years', 'FontSize', 17)


%% Functions

function [dydt] =  SimplifiedEqn(~, Y, nu)

% Load Parameter Values
params = LoadParameters();

% Specify Parameter Values
phiA = params(1);
deltaA = params(2);
lambdaEL = params(3);
omegaEL = params(4);
phiE = params(5);
deltaE = params(6);
C = params(9);
kappa = params(15);

% Relable to easily keep track of compartments
EL = Y(1);
AP = Y(2);
B = Y(3);

dELdt = lambdaEL*AP - omegaEL*(1 - EL/C) - phiE*EL - deltaE*EL;
dAPdt = nu*EL*B - phiA*AP - deltaA*AP;
dBdt = -kappa*EL*B;

dydt = [dELdt; dAPdt; dBdt];

end

%% Event Function
function [value, isterminal, direction] = events(~, y, dydt, IC)
    % B = y(3);
    % value = B - 0.2*IC;
    % disp(['B: ', num2str(B), ', Threshold: ', num2str(0.2 * B_initcond)]);
    value = norm(dydt) - 1e-6;  % Condition: stop when dydt is close to zero
    isterminal = 1;  % Stop the integration
    direction = -1;   % Detect all zeros (both rising and falling)
end